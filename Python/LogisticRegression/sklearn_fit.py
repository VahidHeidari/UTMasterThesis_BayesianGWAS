import matplotlib.pyplot as plt
import numpy as np
import sklearn.linear_model

import fit



def main():
	print('Read data')
	data = fit.ReadData('data.txt')
	x = [ d[0: len(data[0]) - 1] for d in data ]
	y = [ d[len(data[0]) - 1] for d in data ]

	print('Fit sklearn . . .')
	clf = sklearn.linear_model.LogisticRegression().fit(x, y)

	print('Stats:')
	#print('clf      ', clf)
	print('score    ', clf.score(x, y))
	#print('params   ', clf.get_params())
	print('coef     ', clf.coef_)
	print('intercept', clf.intercept_)
	print('classes  ', clf.classes_)
	print('iters    ', clf.n_iter_)

	params = np.concatenate((clf.coef_[0], clf.intercept_))
	print('params :', params)

	tot_true = 0.0
	for i in range(len(data)):
		sig = fit.GetSigmoid(params, data[i][:len(data[0]) - 1])
		y_pred = 1 if sig > 0.5 else 0
		y = data[i][len(data[0]) - 1]
		if y_pred == y:
			tot_true += 1.0
	acc = round(tot_true / len(data) * 100.0, 2)
	print('    Acc:', acc)

	x = [ float(i) / 100.0 for i in range(-8 * 100, 28 * 100) ]
	y = [ fit.GetSigmoid(params, [x[i]]) for i in range(len(x)) ]
	plt.plot(x, y)
	x = [ d[0] for d in data ]
	y = [ d[1] for d in data ]
	plt.scatter(x, y)
	plt.show()


if __name__ == '__main__':
	main()

