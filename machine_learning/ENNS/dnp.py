import copy
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score
import torch
from torch import nn
from typing import List
from warnings import warn


class NeuralNet(nn.Module):
    """
    neural network class, with nn api
    """
    def __init__(self, input_size: int, hidden_size: List[int], output_size: int, regression: bool = False):
        """
        initialization function
        @param input_size: input data dimension
        @param hidden_size: list of hidden layer sizes, arbitrary length
        @param output_size: output data dimension
        @param regression: if true, y is treated as regression response, else classification
        """
        super().__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = output_size
        self.relu = nn.ReLU()
        self.softmax = nn.Softmax(dim=1)
        """layers"""
        self.input = nn.Linear(self.input_size, self.hidden_size[0], bias=False)
        self.hiddens = nn.ModuleList([
            nn.Linear(self.hidden_size[h], self.hidden_size[h + 1]) for h in range(len(self.hidden_size) - 1)])
        self.output = nn.Linear(hidden_size[-1], output_size)
        self.regression = regression

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        forward propagation process, required by the nn.Module class
        @param x: the input data
        @return: the output from neural network
        """
        x = self.input(x)
        x = self.relu(x)
        for hidden in self.hiddens:
            x = hidden(x)
            x = self.relu(x)
        x = self.output(x)
        if not self.regression:
            x = self.softmax(x)
        return x


class DeepNet:
    """
    implements the deep neural network in "https://www.ijcai.org/proceedings/2017/0318.pdf"
    """
    def __init__(self, max_feature: int, num_classes: int = 2, hidden_size: list = None,
                 q: int = 2, num_dropout: int = 50, regression: bool = False, pre_feature: list = []):
        """
        initialization function
        @param max_feature: max number of features selected
        @param num_classes: number of classes in label, ignored if regression is true
        @param hidden_size: a list of hidden layer sizes
        @param q: the norm used in feature selection
        @param num_dropout: number of repetitions of dropout implements when computing gradients
        @param regression: true for regression and false for classification
        @param pre_feature: feature to be included
        """
        self.device = torch.device('cpu')
        if torch.cuda.is_available():
            self.device = torch.device('cuda')
        """model parameters"""
        self.max_feature = max_feature
        self.n = None
        self.p = None
        self.q = q  # norm used in feature selection
        if regression:
            self.num_classes = 1
        else:
            self.num_classes = num_classes
        if hidden_size is None:
            self.hidden_size = [50]
        else:
            self.hidden_size = hidden_size
        self.num_dropout = num_dropout
        self.regression = regression
        """candidate sets"""
        self.pre_feature = pre_feature
        self.S = None
        self.C = None
        self.new = None  # new feature to be added
        self.selected = None
        """model"""
        self.nnmodel = None
        self.last_model = None  # last model
        """optimization parameter"""
        self.criterion = None
        self.learning_rate = None
        self.batch_size = None
        self.epochs = None
        self.dropout_prop = None

    @staticmethod
    def add_bias(x: torch.Tensor) -> torch.Tensor:
        """
        adding bias to input data if the first column is not all ones
        @param x: the input data
        @return: the input data with dummy variable added
        """
        if not all(x[:, 0] == 1):
            x = torch.cat((torch.ones(x.shape[0], 1), x.float()), dim=1)
        return x

    def initialize(self):
        """
        initialize parameters
        """
        if self.pre_feature == []:
            self.S = [0]
            self.selected = []
            self.C = list(range(1, self.p))
        else:
            if 0 not in self.pre_feature:
                self.pre_feature = [0] + self.pre_feature
            self.S = [i for i in self.pre_feature]
            self.selected = [i for i in self.pre_feature if i != 0]
            self.C = [i for i in range(1, self.p) if i not in self.selected]

    def set_parameters(self, x: torch.Tensor, batch_size: int, learning_rate: float, epochs: int,
                       dropout_prop: list = None):
        """
        set optimization parameters in class
        @param x: the input data
        @param batch_size: batch size
        @param learning_rate: learning rate
        @param epochs: number of epochs
        @param dropout_prop: a list of proportions of dropout in each hidden layer
        """
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.epochs = epochs
        if not self.regression:
            self.criterion = nn.CrossEntropyLoss()
        else:
            self.criterion = nn.MSELoss()
        self.n, self.p = x.shape
        if dropout_prop is None:
            self.dropout_prop = [0.5] * len(self.hidden_size)
        else:
            self.dropout_prop = dropout_prop

    def update_nn_weight(self, x: torch.Tensor, y: torch.Tensor, verbosity: int = 0):
        """
        updates the weights of neural network with the selected inputs
        @param x: the training data
        @param y: the training label
        @param verbosity: 0-print everything, >0-print nothing
        """
        input_size = len(self.S)
        hidden_size = self.hidden_size
        output_size = self.num_classes
        self.nnmodel = NeuralNet(input_size, hidden_size, output_size, self.regression).cuda(device=self.device)
        self.xavier_initialization()
        x = x[:, self.S].cuda(device=self.device)
        optimizer = torch.optim.Adam(self.nnmodel.parameters(), lr=self.learning_rate)
        trainset = []
        for i in range(x.shape[0]):
            trainset.append([x[i, :], y[i]])
        trainloader = torch.utils.data.DataLoader(trainset, batch_size=self.batch_size, shuffle=True)
        for e in range(self.epochs):
            running_loss = 0
            for data, label in trainloader:
                input_0 = data.view(data.shape[0], -1).cuda(device=self.device)
                optimizer.zero_grad()
                output = self.nnmodel(input_0.float()).cuda(device=self.device)
                if len(label.shape) > 1:
                    label = label.squeeze(1)
                label = label.cuda(device=self.device)
                if self.regression:
                    output = output.squeeze(1).float().cuda(device=self.device)
                    label = label.float().cuda(device=self.device)
                loss = self.criterion(output, label).cuda(device=self.device)
                loss.backward()
                optimizer.step()
                running_loss += loss.item()
            if e % 10 == 0 and verbosity == 0:
                print(f"Epoch: {e + 1}\nTraining loss: {running_loss/len(trainloader)}")
        self.last_model = copy.deepcopy(self.nnmodel)

    def dropout(self) -> nn.Module:
        """
        returns a copy of randomly dropout model
        """
        if len([p for p in self.dropout_prop if p >= 1 or p <= 0]) > 0:
            raise ValueError("Dropout proportion must be between 0 and 1.")
        if len(self.dropout_prop) > len(self.hidden_size):
            warn(f"Too many dropout proportions, only the first {len(self.hidden_size)} will be used")
        elif len(self.dropout_prop) < len(self.hidden_size):
            warn(f"Too few dropout proportions, dropout won't be applied to the last "
                 f"{len(self.hidden_size) - len(self.dropout_prop) - 1} layers")
        model_dp = copy.deepcopy(self.nnmodel).cuda(device=self.device)
        for h in range(min(len(self.hidden_size) - 1, len(self.dropout_prop))):
            prop = self.dropout_prop[h]
            h_size = self.hidden_size[h]
            dropout_index = np.random.choice(range(h_size), int(h_size * prop), replace=False)
            model_dp.hiddens[h].weight[:, dropout_index] = torch.zeros(
                model_dp.hiddens[h].weight[:, dropout_index].shape).cuda(device=self.device)
        if len(self.hidden_size) <= len(self.dropout_prop):
            prop = self.dropout_prop[len(self.hidden_size) - 1]
            h_size = self.hidden_size[-1]
            dropout_index = np.random.choice(range(h_size), int(h_size * prop), replace=False)
            model_dp.output.weight[:, dropout_index] = torch.zeros(model_dp.output.weight[:, dropout_index].shape)\
                .cuda(device=self.device)
        return model_dp

    def compute_input_gradient(self, x: torch.Tensor, y: torch.Tensor, model: nn.Module) -> torch.Tensor:
        """
        computes the input gradients given a model after dropout
        @param x: the training data
        @param y: the training label
        @param model: the model after dropout
        @return: the input weight gradients
        """
        model_gr = NeuralNet(self.p, self.hidden_size, self.num_classes, self.regression).cuda(device=self.device)
        temp = torch.zeros(model_gr.input.weight.shape).cuda(device=self.device)
        temp[:, self.S] = model.input.weight.cuda(device=self.device)
        model_gr.input.weight.data = temp.cuda(device=self.device)
        for h in range(len(self.hidden_size) - 1):
            model_gr.hiddens[h].weight.data = model.hiddens[h].weight.cuda(device=self.device)
        model_gr.output.weight.data = model.output.weight.cuda(device=self.device)
        output_gr = model_gr(x.float()).cuda(device=self.device)
        if len(y.shape) > 1:
            y = y.squeeze(1)
        y = y.cuda(device=self.device)
        if self.regression:
            output_gr = output_gr.squeeze(1).float().cuda(device=self.device)
            y = y.float().cuda(device=self.device)
        loss_gr = self.criterion(output_gr, y).cuda(device=self.device)
        loss_gr.backward()
        input_gradient = model_gr.input.weight.grad.cuda(device=self.device)
        return input_gradient

    def average_input_gradient(self, x: torch.Tensor, y: torch.Tensor, num_average: int):
        """
        compute the average input gradient over different dropouts
        @param x: the training data
        @param y: the training label
        @param num_average: number of repetitions of dropouts
        @return: the average input weight gradient
        """
        grad_cache = None
        for num in range(num_average):
            model = self.dropout()
            input_grad = self.compute_input_gradient(x, y, model)
            if grad_cache is None:
                grad_cache = input_grad
            else:
                grad_cache += input_grad
        return grad_cache / num_average

    def find_next_input(self, input_gradient: torch.Tensor) -> int:
        """
        computes which input is to be selected next by finding the maximum gradient norm
        @param input_gradient: the average input gradient
        @return: the index of input feature to be included in the model next
        """
        gradient_norm = input_gradient.norm(p=self.q, dim=0)
        gradient_norm = gradient_norm * 1000
        gradient_norm[self.S] = 0
        max_index = torch.argmax(gradient_norm)
        return max_index.item()

    def update_sets(self, max_index):
        """
        updates the selected set and candidate set
        @param max_index: the index of input feature to be included in the model next
        """
        self.S.append(max_index)
        self.C.remove(max_index)
        self.selected = self.S.copy()
        self.selected.remove(0)

    @staticmethod
    def numpy_to_torch(x):
        """
        convert array or dataframe to torch tensor
        """
        if isinstance(x, np.ndarray):
            x = torch.from_numpy(x)
        elif isinstance(x, pd.DataFrame):
            x = torch.from_numpy(x.values)
        return x

    def xavier_initialization(self):
        """
        use Xavier initialization to initialize the newly added feature weight
        @return: None
        """
        if self.last_model is not None:
            weight = torch.zeros(self.hidden_size[0], len(self.S)).cuda(device=self.device)
            nn.init.xavier_normal_(weight, gain=nn.init.calculate_gain('relu'))
            old_s = self.S.copy()
            if self.new in old_s:
                old_s.remove(self.new)
            for i in self.S:
                if i != self.new:
                    weight[:, self.S.index(i)] = self.last_model.input.weight.data[:, old_s.index(i)]\
                        .cuda(device=self.device)
            self.nnmodel.input.weight.data = weight.cuda(device=self.device)
            for h in range(len(self.hidden_size) - 1):
                self.nnmodel.hiddens[h].weight.data = self.last_model.hiddens[h].weight.data.cuda(device=self.device)
            self.nnmodel.output.weight.data = self.last_model.output.weight.data.cuda(device=self.device)
        else:
            print(f"First iter, no Xavier initialization needed.")

    def train(self, x: torch.Tensor, y: torch.Tensor, batch_size: int = 100, learning_rate: float = 0.005,
              epochs: int = 50, dropout_prop: list = None, val: list = None, verbosity: int = 0,
              return_select: bool = False) -> List:
        """
        train the deep neural network on x and y
        @param x: the training data
        @param y: the training label
        @param batch_size: batch size
        @param learning_rate: learning rate
        @param epochs: number of epochs
        @param dropout_prop: a list of proportions of dropout in each hidden layer
        @param val: a list of x_val and y_val
        @param verbosity: 0-print everything; 1-print result only; 2-print nothing
        @param return_select: whether to return the selected variable indices
        """
        x = self.numpy_to_torch(x)
        y = self.numpy_to_torch(y).cuda(device=self.device)
        x = self.add_bias(x).cuda(device=self.device)
        """set parameters"""
        self.set_parameters(x, batch_size, learning_rate, epochs, dropout_prop)
        """initialization"""
        self.initialize()
        """start feature selection"""
        selection = []
        while len(self.S) < self.max_feature + 1:
            self.update_nn_weight(x, y, verbosity=verbosity)
            input_gradient = self.average_input_gradient(x, y, self.num_dropout).cuda(device=self.device)
            self.new = self.find_next_input(input_gradient)
            selection.append(self.new)
            self.update_sets(self.new)
            if verbosity == 0:
                print(f"Number of features selected is {len(self.S) - 1}.")
            if val is not None and verbosity <= 1:
                pred = self.nnmodel(val[0])
                print(f"Validation accuracy is {np.mean(abs(pred.squeeze() - val[1].squeeze()))}")
        self.update_nn_weight(x, y, verbosity=verbosity)
        if verbosity <= 1:
            print(f"Feature selection completed, selected features are {self.selected}")
        if return_select:
            return selection
        else:
            return []

    def train_return_next(self, x: torch.Tensor, y: torch.Tensor, batch_size: int = 100, learning_rate: float = 0.005,
                         epochs: int = 50, dropout_prop: list = None, val: list = None, verbosity: int = 0,
                         return_select: bool = False, initial: list = None) -> List:
        """
        train with some initially added predictors and return the next
        @param x: the training data
        @param y: the training label
        @param batch_size: batch size
        @param learning_rate: learning rate
        @param epochs: number of epochs
        @param dropout_prop: a list of proportions of dropout in each hidden layer
        @param val: a list of x_val and y_val
        @param verbosity: 0-print everything; 1-print result only; 2-print nothing
        @param return_select: whether to return the selected variable indices
        @param initial: the initially added predictors
        """
        if initial is None:
            self.train(x, y, batch_size, learning_rate, epochs, dropout_prop, val, verbosity, return_select)
        else:
            x = self.numpy_to_torch(x)
            y = self.numpy_to_torch(y)
            x = self.add_bias(x)
            self.S = initial
            if 0 not in self.S:
                self.S = [0] + self.S
            input_size = len(self.S)
            hidden_size = self.hidden_size
            output_size = self.num_classes
            self.nnmodel = NeuralNet(input_size, hidden_size, output_size, self.regression)
            self.set_parameters(x, batch_size, learning_rate, epochs, dropout_prop)
            self.update_nn_weight(x, y, verbosity)
            self.last_model = copy.deepcopy(self.nnmodel)
            self.update_nn_weight(x, y, verbosity)
            input_gradient = self.average_input_gradient(x, y, self.num_dropout)
            self.new = self.find_next_input(input_gradient)
            return [self.new]


    def predict(self, x, prob: bool = False) -> torch.Tensor:
        """
        making prediction with the trained model of the given x
        @param x: the testing data
        @param prob: whether to return probability or class labels, ignored if regression is true
        @return: the prediction of y
        """
        if self.nnmodel is None:
            raise ValueError("Model not trained, please run train first.")
        x = self.numpy_to_torch(x)
        x = self.add_bias(x).cuda(device=self.device)
        if x.shape[1] != self.p:
            raise ValueError("Dimension of x is wrong.")
        x = x[:, list(self.S)].float().cuda(device=self.device)
        y_pred = self.nnmodel(x)
        if self.regression or prob:
            return y_pred
        y_pred = torch.argmax(y_pred, dim=1)
        return y_pred

    def predict_error(self, x, y, return_res: bool = False) -> (float, float, float):
        """
        makes prediction on x and computes the prediction error from y
        @param x: the testing data
        @param y: the testing label
        @param return_res: whether to return results or just print them
        @return: the accuracy, auc and f1_score
        """
        prob_pred = self.predict(x, prob=True)
        prob_pred = prob_pred[:, 1].squeeze().float()
        y_pred = self.predict(x)
        y_pred = y_pred.squeeze().float()
        y = self.numpy_to_torch(y)
        y = y.squeeze().float()
        if not self.regression:
            acc = 1 - (abs(y - y_pred)).mean().item()
            auc = roc_auc_score(y.detach().numpy(), prob_pred.detach().numpy())
            f1 = f1_score(y.detach().numpy(), y_pred.detach().numpy())
            print(f"Testing accuracy is {acc}.\nTesting auc is {auc}.\nTesting f1 score is {f1}.")
        else:
            mse = np.mean(np.square(y_pred.numpy() - y.numpy()))
            mae = np.mean(np.abs(y_pred.numpy() - y.numpy()))
            mape = np.mean(np.abs(y_pred.numpy() - y.numpy()) / np.abs(y.numpy()))
        if return_res:
            if not self.regression:
                return acc, auc, f1
            else:
                return mse, mae, mape
