from LastProject.simulation import generate_data, linear_data
from LastProject.TwoStepNN import TwoStepNet
import numpy as np
from LastProject.dnp import DeepNet
from sklearn.linear_model import LogisticRegression, Lasso
from pyHSICLasso import HSICLasso



def generateXY(seed: int, struct: str = 'Linear', response: str = 'classification', n: int = 1000, p: int = 10000,
               s: int = 5, mu: float = 0, sigma: float = 0, error: float = 1) -> (np.array, np.array):
    """generates data"""
    if struct == 'NeuralNet':
        if response == 'classification':
            x, _, y, _ = generate_data(seed, s, n * 2, n, 0, mu, sigma, False)
        elif response == 'regression':
            x, _, y, _ = generate_data(seed, s, n * 2, n, 0, mu, sigma, True)
            y = y.ravel() + np.random.randn(n) * error
        else:
            print("response need to be classification or regression")
            return None, None
    elif struct in ['Linear', 'Additive']:
        x = np.random.rand(n, p).reshape(n, p) * 2 - 1
        if struct == 'Linear':
            beta = np.random.randn(s)
            beta[abs(beta) < 1e-1] = 1e-1
            eta = np.matmul(x[:, :s], beta)
        else:
            eta = np.sin(x[:, 0]) + x[:, 1] + np.exp(x[:, 2]) + np.square(x[:, 3]) + np.log(x[:, 4] + 2) - 2
        if response == 'classification':
            prob = 1 / (1 + np.exp(-eta))
            y = np.array([np.random.binomial(1, prob[i]) for i in range(n)]).reshape(-1, 1)
        elif response == 'regression':
            y = eta + np.random.randn(n) * error
        else:
            print("response need to be classification or regression")
            return None, None
    return x, y




def right_select(true, select):
    count_tp = 0
    count_fp = 0
    if 0 in select:
        select.remove(0)
    for i in range(len(select)):
        if select[i] in true:
            count_tp += 1
        else:
            count_fp += 1
    return count_tp, count_fp


def lasso_selection(x, y, num: int, response: str = "Classification") -> list:
    c1 = 1
    c2 = 1
    y = y.ravel()
    if response == "Classification":
        model1 = LogisticRegression(penalty='l1', solver='liblinear', C=c1)
        model1.fit(x, y)
        coefs = model1.coef_.ravel()
        number_of_selection = len([i for i in coefs if i != 0])
        if number_of_selection == num:
            return sorted([i+1 for i, j in enumerate(coefs) if j != 0])
        while number_of_selection > num:
            print(f"Number of selected variables is {number_of_selection}. {c1}")
            c1 /= 2
            model1 = LogisticRegression(penalty='l1', solver='liblinear', C=c1)
            model1.fit(x, y)
            coefs = model1.coef_.ravel()
            number_of_selection = len([i for i in coefs if i != 0])
        model1 = LogisticRegression(penalty='l1', solver='liblinear', C=c2)
        model1.fit(x, y)
        coefs = model1.coef_.ravel()
        number_of_selection = len([i for i in coefs if i != 0])
        while number_of_selection < num:
            print(f"Number of selected variables is {number_of_selection}.")
            c2 *= 2
            model1 = LogisticRegression(penalty='l1', solver='liblinear', C=c2)
            model1.fit(x, y)
            coefs = model1.coef_.ravel()
            number_of_selection = len([i for i in coefs if i != 0])
        while number_of_selection != num:
            print(f"Number of selected variables is {number_of_selection}.")
            c = (c1 + c2) / 2
            model1 = LogisticRegression(penalty='l1', solver='liblinear', C=c)
            model1.fit(x, y)
            coefs = model1.coef_.ravel()
            number_of_selection = len([i for i in coefs if i != 0])
            if number_of_selection > num:
                c2 = c
            else:
                c1 = c
        return sorted([i+1 for i, j in enumerate(coefs) if j != 0])
    else:
        model1 = Lasso(alpha=c1)
        model1.fit(x, y)
        coefs = model1.coef_.ravel()
        number_of_selection = len([i for i in coefs if i != 0])
        if number_of_selection == num:
            return sorted([i+1 for i, j in enumerate(coefs) if j != 0])
        while number_of_selection > num:
            print(f"Number of selected variables is {number_of_selection}. {c1}")
            c1 *= 2
            model1 = Lasso(alpha=c1)
            model1.fit(x, y)
            coefs = model1.coef_.ravel()
            number_of_selection = len([i for i in coefs if i != 0])
        model1 = Lasso(alpha=c2)
        model1.fit(x, y)
        coefs = model1.coef_.ravel()
        number_of_selection = len([i for i in coefs if i != 0])
        while number_of_selection < num:
            print(f"Number of selected variables is {number_of_selection}.")
            c2 /= 2
            model1 = Lasso(alpha=c2)
            model1.fit(x, y)
            coefs = model1.coef_.ravel()
            number_of_selection = len([i for i in coefs if i != 0])
        while number_of_selection != num:
            print(f"Number of selected variables is {number_of_selection}.")
            c = (c1 + c2) / 2
            model1 = Lasso(alpha=c)
            model1.fit(x, y)
            coefs = model1.coef_.ravel()
            number_of_selection = len([i for i in coefs if i != 0])
            if number_of_selection > num:
                c2 = c
            else:
                c1 = c
        return sorted([i+1 for i, j in enumerate(coefs) if j != 0])




