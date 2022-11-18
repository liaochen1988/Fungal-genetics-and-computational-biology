import numpy as np


def relu(x):
	return np.maximum(x, 0)


def sigmoid(x):
	return np.exp(x) / (1 + np.exp(x))


def generate_data(
		seed: int, m: int, n: int = 8500, train_size: int = 1000, flip: float = 0.05, mu: float = 0, sig: float = 1,
		regression: bool = False) -> (np.array, np.array, np.array, np.array):
	"""
	generate random simulation data according to the set up in "https://www.ijcai.org/proceedings/2017/0318.pdf"
	@param train_size: training sample size
	@param n: total sample size
	@param m: nonzero input
	@param seed: seed
	@param flip: proportion of labels to be flipped
	@param mu: the mean of signal
	@param sig: the variance of nonzero weights
	@param regression: true for regression and false for classification
	@return: x_train, x_test, y_train, y_test
	"""
	"""parameters"""
	np.random.seed(seed)
	p = 10000
	h1 = 50
	h2 = 30
	h3 = 15
	h4 = 10
	sig1 = np.sqrt(1)
	sig = np.sqrt(sig)
	"""generate x"""
	x = np.random.rand(n*p).reshape(n, p) * 2 - 1
	"""neural network forward"""
	w1 = np.concatenate([np.random.randn(m, h1) * sig + mu, np.zeros([p-m, h1])])
	w2 = np.random.randn(h1, h2) * sig1
	w3 = np.random.randn(h2, h3) * sig1
	w4 = np.random.randn(h3, h4) * sig1
	w5 = np.random.randn(h4, 1) * sig1
	z1 = np.matmul(x, w1)
	a1 = relu(z1)
	z2 = np.matmul(a1, w2)
	a2 = relu(z2)
	z3 = np.matmul(a2, w3)
	a3 = relu(z3)
	z4 = np.matmul(a3, w4)
	a4 = relu(z4)
	if not regression:
		prob = sigmoid(np.matmul(a4, w5))
		"""generate y"""
		y = np.where(prob > 0.5, 1, 0)
		"""flip labels"""
		flip = np.random.choice(range(n), size=int(n*flip))
		y[flip, :] = 1 - y[flip, :]
		print(f"y has mean {y.mean()}")
	else:
		y = np.matmul(a4, w5)
		y = y.reshape(-1) + np.random.randn(n)
	return x[:train_size, :], x[train_size:, :], y[:train_size], y[train_size:]


def linear_data(
		seed: int, m: int, n: int = 8500, train_size: int = 1000, flip: float = 0.05, mu: float = 0, sig: float = 1,
		lb: float = 0.5) -> (np.array, np.array, np.array, np.array):
	"""
	@param train_size: training sample size
	@param n: total sample size
	@param m: nonzero input
	@param seed: seed
	@param flip: proportion of labels to be flipped
	@param mu: the mean of signal
	@param sig: the variance of nonzero weights
	@param lb: lower bound of signal strength
	@return: x_train, x_test, y_train, y_test
	"""
	"""parameters"""
	np.random.seed(seed)
	p = 10000
	sig = np.sqrt(sig)
	"""generate x"""
	x = np.random.rand(n*p).reshape(n, p) * 2 - 1
	"""generate beta"""
	lb = abs(lb)
	beta = np.concatenate([np.random.randn(m) * sig + mu, np.zeros(p - m)]).reshape(p, 1)
	beta[(0 < beta) & (beta < lb * sig)] = lb * sig
	beta[(-lb * sig < beta) & (beta < 0)] = -lb * sig
	print(f"beta is {beta[:m]}")
	"""calculate probability"""
	eta = np.matmul(x, beta)
	prob = sigmoid(eta)
	"""generate y"""
	y = np.array([np.random.binomial(1, prob[i]) for i in range(n)]).reshape(-1, 1)
	"""flip labels"""
	flip = np.random.choice(range(n), size=int(n*flip))
	y[flip, :] = 1 - y[flip, :]
	print(f"y has mean {y.mean()}")
	return x[:train_size, :], x[train_size:, :], y[:train_size], y[train_size:]