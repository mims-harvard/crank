import numpy as np
from scipy import stats


def _smooth(x):
    smooth = x.copy()
    n_changes = 1
    while n_changes != 0:
        prev = smooth.copy()
        for i in xrange(1, len(x)-1):
            smooth[i] = np.median(prev[i-1:i+2])
        smooth[0] = np.median([prev[0], smooth[1], 3 * smooth[1] - 2 * smooth[2]])
        smooth[-1] = np.median([prev[-1], smooth[-2], 3 * smooth[-2] - 2 * smooth[-3]])
        n_changes = np.sum(smooth != prev)
    return smooth


def _cummax(x):
    y = np.array([np.max(x[:i]) for i in xrange(1, len(x)+1)])
    return y


def rra(data, prior=0.2, num_bin=4, num_iter=10,
        return_all=True, corr_stop=1):
    # data:
    #   - data is a numpy matrix with objects in rows and different ranked lists in columns.
    #   - data is a real-valued matrix where entries with higher values indicate stronger preference.
    #
    # Note 1: Note that a higher value in data matrix is better and indicates a higher-priority of an object. When
    # converted to ranks the largest value gets rank 1. If input data are ranks (i.e., lower values indicate higher
    # priority), the ranks might need to be reversed.
    nr, nc = data.shape
    nrp = int(np.floor(nr * prior))
    print 'Nrp: %d' % nrp
    rank_data = np.zeros(data.shape)
    for i in xrange(nc):
        rank_data[:, i] = stats.rankdata(-data[:, i]) / float(nr)

    bayes_factors = np.zeros((num_bin, nc))
    binned_data = np.ceil(rank_data * num_bin).astype('int')
    bayes_data = np.zeros((nr, nc))

    # estimated ranks, smaller is better (= closer to the top)
    guess = np.mean(rank_data, 1)
    cprev = 0
    for iter in xrange(num_iter):
        if corr_stop - cprev < 1e-15:
            print 'Converged!'
            break

        # assign the top np of aggregated predictions to be the positive class and the
        # rest of the predictions to the negative class
        guess_last = guess.copy()
        oo = np.argsort(guess)
        guess[oo[:nrp]] = 1.
        guess[oo[nrp:]] = 0.

        # Heuristic to make the approach more robust:
        # -- computing Bayes factors cumulatively starting from the top bin
        for i in xrange(nc):
            for bin in xrange(1, num_bin + 1):
                tpr = np.sum(guess[binned_data[:, i] <= bin])
                fpr = np.sum((1. - guess)[binned_data[:, i] <= bin])
                bayes_factors[bin-1, i] = np.log((tpr + 1.) / (fpr + 1.) / (prior / (1. - prior)))

        # Heuristic to make the approach more robust:
        # -- smoothing using Tukey's running median
        # -- enforcing a monotone decrease of Bayes factors
        for i in xrange(nc):
            bayes_factors[:, i] = _smooth(bayes_factors[:, i])
            bayes_factors[:, i] = _cummax(bayes_factors[:, i][::-1])[::-1]

        # winzorization step: for each bin we decrease the maximum Bayes factor
        # to that of the second maximum
        # for bin in xrange(num_bin):
        #     oo = np.argsort(-bayes_factors[bin, :])
        #     bayes_factors[bin, oo[0]] = bayes_factors[bin, oo[1]]

        for i in xrange(nc):
            # if j-th entry in i-th ranking falls into k-th bin,
            # then bayes_data[j, i] should be the bayes factor of k-th
            # bin in i-th ranking
            bayes_data[:, i] = bayes_factors[binned_data[:, i]-1, i]

        # bb = np.exp(np.sum(bayes_data, 1))
        # f = prior / (1. - prior)
        # prob = bb * f / (1. + bb * f)
        # exp = np.sort(-prob)[nrp]

        guess = stats.rankdata(-np.sum(bayes_data, 1))
        cprev = stats.pearsonr(guess, guess_last)[0]
        print 'Correlation with previous iteration: %f' % cprev
    if return_all:
        return guess, bayes_data, bayes_factors
    else:
        return guess
