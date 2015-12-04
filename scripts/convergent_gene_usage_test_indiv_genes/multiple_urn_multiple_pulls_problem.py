import sys
import itertools

def prob_with_multiple_pulls(prob_of_sel, pulls):
    """This subscript is the part that models the multiple pulls drawn from each urn (with replacement). 'prob_of_sel' = probability of selecting a brown turd; 'pulls' = number of pulls made from the urn.
    Specifically, this models the probability of pulling a brown turd from an urn at least once."""
    prob = 1
    for i in xrange(pulls):
        prob *= 1 - prob_of_sel
    return 1 - prob

def run(list_o_probs, num_of_pulls, only_max_num_urns=None):
    """This script gives the pdf for a multile urn problem. Specifically, each urn has a probability of selecting a browm turd. These probs for each urn are listed in 'list_o_probs'. If a selction is made from each urn, this script provides the probability that X number of brown turds will be drawn, where 0<=X<=[number of urns]."""
    """This script has been adapted from 'multiple_urn_problem.py', in order to additionally model multiple pulls being drawn from each of the urns (this is a mre realistic model for the influenza vaccination stuff). So each urn has a number of times that independent draws are made from it (with replacement) ('num_of_pulls)."""
    """UPDATE: We added the parameter 'only_max_num_urns', and when this is defined (could be anything) then this script only returns the probability of pulling a brown turd from each urn. That is the maximum number of urns return a brown turd. We did this because when running the power calculations using 'power_across_many_parameters_gene_by_gene_convergence_test.py', in this case returnig the entire pdf was taking too long, and we were only interested in the case where all urns yeild a turd (i.e. gene is convergent) so we changed this script so that it only returns this value."""
    pdf = []
    index_of_urns = range(len(list_o_probs))
    #for each possible number of brown turds
    for i in xrange(len(list_o_probs)+1):
        #if we are to only return the probabiility of pulling a brown
        #turd from the max number of urns and 'i' is not at the max 
        #number of urns, then skip
        if only_max_num_urns and i != len(list_o_probs):
            continue
        if i == 0:
            prob = 1
            for j in xrange(len(list_o_probs)):
                prob *= 1 - prob_with_multiple_pulls(list_o_probs[j], num_of_pulls[j])
        elif i == len(list_o_probs):
            prob = 1
            for j in xrange(len(list_o_probs)):
                prob *= prob_with_multiple_pulls(list_o_probs[j], num_of_pulls[j])
            #if we are to to only return the probabiility of pulling a brown
            #for the max number of urns, then just return this prob and quit
            if only_max_num_urns:
                return prob
        else:
            prob = 0
            for j in itertools.combinations(index_of_urns, i):
                prob_subset = 1
                for k in index_of_urns:
                    if k in j:
                        prob_subset *= prob_with_multiple_pulls(list_o_probs[k], num_of_pulls[k])
                    else:
                        prob_subset *= (1 - prob_with_multiple_pulls(list_o_probs[k], num_of_pulls[k]))
                prob += prob_subset
        pdf.append(prob)
    return pdf

if __name__ == '__main__':
    lists = sys.argv[1:]
    list_o_probs = lists[:len(lists)/2]
    num_of_pulls = lists[len(lists)/2:]
    list_o_probs = [float(i) for i in list_o_probs]
    num_of_pulls = [int(i) for i in num_of_pulls]
    run(list_o_probs, num_of_pulls)
