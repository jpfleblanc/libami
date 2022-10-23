import pyami_example as ex
import pyami
import time

def main():
    example_2()

def example_2():

    print("-_-_-_ Example - Second Order _-_-_-")

    # start example
    print("\n-----Constructing AMI SPR format -----\n")
    
    # class instance
    ami = pyami.AmiBase()

    # second order self energy setup
    R0 = ex.construct_example2()

    # external variables
    avars = ex.construct_ext_example2()

    # timing info for setup
    t1 = time.time()

    # These are opaque types so that they are mutable in c++ construct function.
    S_array = pyami.S_t()
    P_array = pyami.P_t()
    R_array = pyami.R_t()

    # Integration/Evaluation parameters
    E_REG = 0 # numberical regulator for small energies.  If inf/nan results try E_REG=1e-8 
    N_INT = 2 # number of matsubara sums to perform 
    test_amiparms = pyami.AmiBase.ami_parms(N_INT, E_REG)

    # now construct!
    ami.construct(test_amiparms, R0, R_array, P_array, S_array)

    # time construct leg and start evaluate time
    t2 = time.time()

    # Evaluate the integrand for ext parms in avars
    calc_result = ami.evaluate(test_amiparms, R_array, P_array, S_array, avars)
    
    # time the end of eval
    t_end = time.time()

    # millisecond result
    diff1 = (t2 - t1) * 1000
    diff2 = (t_end - t2) * 1000

    # print results
    print(f"Result was {calc_result}")
    print(f"Construction took {diff1} microseconds")
    print(f"Evaluation took {diff2} microseconds")

    # end example

    # start problem with optimization
    print("\n-----Constructing AMI Optimized SPR format -----\n")

    t3 = time.time()

    unique = pyami.g_prod_t()
    rref = pyami.R_ref_t()
    eval_list = pyami.R_ref_t()

    #take existing solution from first part and factorize it 
    ami.factorize_Rn(R_array[len(R_array)-1], unique, rref, eval_list)

    # end factor time
    t4 = time.time()

    opt_calc_result = ami.evaluate(test_amiparms, R_array, P_array, S_array, avars, unique, rref, eval_list)

    # end eval time
    t5 = time.time()

    #timing stuff in milliseconds
    diff3 = (t4 - t3) * 1000
    diff4 = (t5 - t4) * 1000

    print(f"Optimized result was {opt_calc_result}")
    print(f"Factorization returned {len(unique)} unique Green's functions and took {diff3} microseconds")
    print(f"Evaluation took {diff4} microseconds")

    # end opt example
    # start terms example
    print("\n-----Constructing AMI term by term-----\n")

    t6 = time.time()

    # simplified storage types
    amiterms = pyami.terms()

    # construct soln for prob in R0
    ami.construct(N_INT, R0, amiterms)

    t7 = time.time()

    #Evaluate the term-by-term result for external values in 'avars'. Note that the test_amiparms is the same as the first case
    term_val = ami.evaluate(test_amiparms, amiterms, avars)

    t8 = time.time()

    diff5 = (t7 - t6) * 1000
    diff6 = (t8 - t7) * 1000
    print(f"Construction took {diff5} microseconds")
    print(f"Result has num_terms={len(amiterms)} compared to standard {len(R_array[N_INT])}")
    print(f"Term result was {term_val}")
    print(f"Evaluation took {diff6} microseconds")

    # end term example
    # start optimized term example
    print("\n-----Constructing AMI Optimized term by term-----\n")

    t9 = time.time()

    ami.factorize_terms(amiterms, unique, rref, eval_list)

    t10 = time.time()

    # Evaluate optimized term-by-term construction 
    opt_term_val = ami.evaluate(test_amiparms, amiterms, avars, unique, rref, eval_list)

    t11 = time.time()
    diff7 = (t9 - t8) * 1000
    diff8 = (t10 - t9) * 1000

    print(f"OPT term val= {opt_term_val}")
    print(f"Factorization returned {len(unique)} unique Green's functions and took {diff7} microseconds")
    print(f"Evaluation took {diff8} microseconds")

if __name__ == '__main__':
    main()


