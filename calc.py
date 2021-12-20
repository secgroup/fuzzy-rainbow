import math

import matplotlib
import numpy as np
import pylab

# USAGE:
# To reproduce the parameters and the estimated performance for DES and A5/1 
# presented in the paper set the following variable to one of the following 
# and execute the tool:

CIPHER = 'DES'    # (Table 1)
# CIPHER = 'A5/1 90%'  # (Table 2, first column)
# CIPHER = 'A5/1 99%'  # (Table 2, second column)

#
# To use the tool for a specific application, define a new set of global 
# parameters with hardware and attack constraints. 
# Tune compression parameters by adjusting the values based on 
# memory constraints and the tool output.
#
# NOTE: the tool automatically finds the optimal s value starting from 
# value 15 and iterating the search until a fixpoint is reached.
#

if CIPHER == 'DES':
    # ============ DES parameters (Table 1) ============
    FREQ = 115 * 1000 * 1000  # FPGA frequency
    CORES = 1664  # number of compute units
    CLOCKS = 4  # number of clock required to perform encryption
    ENCRATE = 2 ** 35.48  # false if parameters above used else only encryption rate is used to compute FPGA time
    # disks
    NDISKS = 8
    IOPSSINGLE = 300000  # IOPS found with FIO, better if cumulative with NDISK, i.e. IOPSSSINGLE = FIOIOPS/NDISKS
    ADJFACTOR = 0.7  # adjustment factor
    # optimization parameters
    FMSCRANGE = (0.3, 5, 200)  # Fmsc ranges for optimization (min,max,numbe of points)
    FMSCMINOFFSET = 0  # offset from min value in Fpc vs Ftc,s; positive right, negative left
    FIGURE = False  # show Fpc vs Ftc,s figure
    FINDBESTS = True  # find best S values
    SRANGE = (15, 101)
    # input params
    NUMPAIRS = 1
    TARGETSR = 0.99  # target success rate
    BRUTEFORCE = 0  # if bruteforce needed in search
    NSPACESIZE = 56  # search space size in bits
    BLOCKSIZE = 64  # block size, in case it differs from log_2(N), e.g. block cipher DES
    CHOSENMEMORY = 40.84  # Memory size (log2), change here to adjust time/memory tradeoff
    # compression params
    TRUNCATIONBITS = 10  # bits to truncate -> calculate r as size of ep - DPbits - truncation bits
    INDEXFILECOMPBITS = 6  # bits to remove with index file method on EP

elif CIPHER == 'A5/1 90%':
    # ============ A5/1 parameters (Table 2, first column) ============
    FREQ = 115 * 1000 * 1000 # FPGA frequency
    CORES = 1664 # number of compute units
    CLOCKS = 4 # number of clock required to perform encryption
    ENCRATE = 2**35.48/(164/64) # false if parameters above used else only encryption rate is used to compute FPGA time
    # disks
    NDISKS = 8
    IOPSSINGLE = 300000 # IOPS found with FIO, better if cumulative with NDISK, i.e. IOPSSSINGLE = FIOIOPS/NDISKS
    ADJFACTOR = 0.7 # adjustment factor
    # optimization parameters
    FMSCRANGE = (0.3,5,200) # Fmsc ranges for optimization (min,max,numbe of points)
    FMSCMINOFFSET=0 # offset from min value in Fpc vs Ftc,s; positive right, negative left
    FIGURE = False # show Fpc vs Ftc,s figure
    FINDBESTS = True # find best S values
    SRANGE = (15,101)
    # input params
    NUMPAIRS = 51
    TARGETSR = 0.90 # target success rate
    BRUTEFORCE = 0 # if bruteforce needed in search
    NSPACESIZE = 61.36 # search space size in bits
    BLOCKSIZE = 64 # block size, in case it differs from log_2(N), e.g. block cipher DES
    CHOSENMEMORY = 40.84  # Memory size (log2), change here to adjust time/memory tradeoff
    # compression params
    TRUNCATIONBITS = 15 # bits to truncate -> calculate r as size of ep - DPbits - truncation bits
    INDEXFILECOMPBITS = 9 # bits to remove with index file method on EP

elif CIPHER == 'A5/1 99%':
    # ============ A5/1 parameters (Table 2, second column) ============
    FREQ = 115 * 1000 * 1000 # FPGA frequency
    CORES = 1664 # number of compute units
    CLOCKS = 4 # number of clock required to perform encryption
    ENCRATE = 2**35.48/(164/64) # false if parameters above used else only encryption rate is used to compute FPGA time
    # disks
    NDISKS = 8
    IOPSSINGLE = 300000 # IOPS found with FIO, better if cumulative with NDISK, i.e. IOPSSSINGLE = FIOIOPS/NDISKS
    ADJFACTOR = 0.7 # adjustment factor
    # optimization parameters
    FMSCRANGE = (0.3,5,200) # Fmsc ranges for optimization (min,max,numbe of points)
    FMSCMINOFFSET=0 # offset from min value in Fpc vs Ftc,s; positive right, negative left
    FIGURE = False # show Fpc vs Ftc,s figure
    FINDBESTS = True # find best S values
    SRANGE = (15,101)
    # input params
    NUMPAIRS = 204
    TARGETSR = 0.99 # target success rate
    BRUTEFORCE = 0 # if bruteforce needed in search
    NSPACESIZE = 61.36 # search space size in bits
    BLOCKSIZE = 64 # block size, in case it differs from log_2(N), e.g. block cipher DES
    CHOSENMEMORY = 40.84  # Memory size (log2), change here to adjust time/memory tradeoff
    # compression params
    TRUNCATIONBITS = 15 # bits to truncate -> calculate r as size of ep - DPbits - truncation bits
    INDEXFILECOMPBITS = 12 # bits to remove with index file method on EP
else:
    print('[-] Unsupported cipher, check the CIPHER variable')
    exit(1)

print('[*] Cipher: {}\n'.format(CIPHER))
# model from Byoung-Il Kim and Jin Hong "Analysis of the Non-Perfect Table Fuzzy Rainbow Tradeoff"
# mtls N and r are globals, factor pairs (number of chal/resp pairs) added to the model to extend it
# to the multi-target version of the algorithm

# matrix stopping constant
def Fmsc():
    return (m * t ** 2 * s) / N


# expected number of i-th color boundary points in a fuzzy rainbow matrix
def m_i(i):
    if i == 0:
        return m
    return 2 * m / (2 + Fmsc() * (i / s))


# coverage rate of a fuzzy rainbow matrix
def Fcr():
    return (2 / Fmsc()) * math.log(1 + Fmsc() / 2)


# coverage rate of a fuzzy rainbow matrix alternative formula
def Fcr2():
    return 1 / (m * t * s) * sum(DM_i(i) for i in range(s))


# expecter number of distinct points in the i-th color DP submatrix
def DM_i(i):
    return m_i(i) * t


# precomputation coefficient
def Fpc():
    return (m * t * s * l) / N


# precomputation coefficient alternative formula
def Fpc_alt_2():
    return (-math.log(1 - Fps())) / Fcr()


# success rate of inversion having multiple data pairs (chal/resp pairs)
def Fps_multi(pairs):
    return 1 - pow(math.exp(-Fcr() * Fpc()), pairs)


# succcess rate of inversion
def Fps():
    return 1 - math.exp(-Fcr() * Fpc())


# tradeoff coefficient for the fuzzy rainbow tradeoff. Part of tradeoff curve TM^2 = F_{tc,s} N^2
def Ftc_s(s):
    return (
        Fmsc() ** 2
        * (l / t) ** 3
        * sum(
            ((1 - (i - 1) / s) * (1 + Fmsc() * (i / s)) + (Fmsc() / (s ** 2)))
            * ((2 + Fmsc() * i / s) / (2 + Fmsc())) ** (2 * l / t)
            * (1 / s)
            for i in range(1, int(s) + 1)
        )
    )


# total number of iterations of f for generation of online chains
def Tf(pairs):
    return (
        pairs
        * t
        * l
        * sum((s - i + 1) * ((2 + Fmsc() * i / s) / (2 + Fmsc())) ** (2 * l / t) for i in range(1, int(s) + 1))
    )


# total number of iterations of f for resolving alarms during online phase
def Ta(pairs):
    return (
        pairs
        * (t * l * Fmsc() / s)
        * sum(
            (i * (s - i + 1) + 1) * (((2 + Fmsc() * i / s) / (2 + Fmsc())) ** (2 * l / t))
            for i in range(1, int(s) + 1)
        )
    )


# expected number of table lookups for the online phase
def Lookups(pairs):
    return pairs * l * sum((((2 + Fmsc() * i / s) / (2 + Fmsc())) ** (2 * l / t)) for i in range(1, int(s) + 1))


# extra invocations of f for resolving truncation-related alarms during online phase
def Ftrunc(pairs):
    return (
        pairs
        * t
        * l
        * (m / r)
        * sum(i * ((2 + Fmsc() * i / s) / (2 + Fmsc())) ** (2 * t / l) for i in range(1, int(s) + 1))
    )


# upper bound number of table lookups for the online phase
def LookupsUpperBound(pairs):
    return pairs * l * s


# Wall clock time formulas

# FPGA computation time
def FPGATime(F):
    if ENCRATE is not False:
        return F / ENCRATE
    return F / FREQ / CORES * CLOCKS
    # return (F * CLOCKS)/(FREQ * CORES)


# FPGA encryption rate
def CalcFPGAEncRate():
    if ENCRATE is not False:
        return ENCRATE
    return (FREQ * CORES) / CLOCKS


# Disk access time
def DiskTime(Lookups):
    systemIOPS = NDISKS * IOPSSINGLE * ADJFACTOR
    return Lookups / systemIOPS


# formulas to compute optimal values from Fps and s

# compute success rate for finding optimal values with kim and hong method starting
# from target success rate and number of data pairs
def Fps_from_multi(pairs, Fpsmulti):
    return 1 - math.pow((1 - Fpsmulti), (1 / pairs))


# ratio between l and t parameters of the tradeoff as function of variable Fmsc
def ltratio(Fmsc, Fps):
    # return Fpc/Fmsc
    return (-math.log(1 - Fps)) / (Fmsc * Fcr_fixedmsc(Fmsc))


# tradeoff coefficient for the fuzzy rainbow tradeoff. Part of tradeoff curve TM^2 = F_{tc,s} N^2
# formula using ratio between l and t, which can be seen as a function of single variable Fmsc
def Ftc_s_ratio(s, ltratio, Fmsc):
    return (
        Fmsc ** 2
        * (ltratio) ** 3
        * sum(
            ((1 - (i - 1) / s) * (1 + Fmsc * (i / s)) + (Fmsc / (s ** 2)))
            * ((2 + Fmsc * i / s) / (2 + Fmsc)) ** (2 * ltratio)
            * (1 / s)
            for i in range(1, int(s) + 1)
        )
    )


# coverage rate of a fuzzy rainbow matrix as function of single variable Fmsc
def Fcr_fixedmsc(Fmsc):  # coverage rate
    return (2 / Fmsc) * math.log(1 + Fmsc / 2)


# precomputation coefficient as function of single variable Fmsc
def Fpc_alt(Fps, Fmsc):
    return (-math.log(1 - Fps)) / Fcr_fixedmsc(Fmsc)


# formulas to compute m,t,l parameters once a point in the fpc vs ftcs curve is chosen and a (T,M)-pair
# satisfying the tradeoff curve is selected


def smallm(Fmsc, Fcr, Ftcs, s, Fps, N, T):
    return (Fmsc * math.pow(Fcr, 2) * Ftcs * s * N) / (math.pow(math.log(1 - Fps), 2) * T)


def smallt(Fps, Fcr, Ftcs, s, T):
    return -(math.log(1 - Fps)) / (Fcr * math.sqrt(Ftcs) * s) * math.sqrt(T)


def smalll(Fps, Fmsc, Fcr, Ftcs, s, T):
    return math.pow(math.log(1 - Fps), 2) / (Fmsc * math.pow(Fcr, 2) * math.sqrt(Ftcs) * s) * math.sqrt(T)


smin = 16
SVALUE = 1000
round = 0
# Searches for the optimal value of s
while smin != SVALUE:
    SVALUE = smin
    round+=1
    print("======================================================")
    print("Iteration {}, s = {}".format(round,SVALUE))
    print("======================================================")
    # Finding optimal parameters from input success rate, number of pairs.
    print("[-] Finding optimal parameters...")
    Fps_fix = Fps_from_multi(NUMPAIRS, TARGETSR)  # target success rate
    s = SVALUE  # between 15 and 100 #

    Fmscrange = np.linspace(*FMSCRANGE)

    Fpcrange = [Fpc_alt(Fps_fix, Fmscel) for Fmscel in Fmscrange]
    Ftcsrange = [Ftc_s_ratio(s, ltratio(Fmscel, Fps_fix), Fmscel) for Fmscel in Fmscrange]

    xmin = np.argmin(Ftcsrange)

    # print(
    #     "    xmin: {:d}, Fpcmin: {:.3f}, Ftcsmin: {:.3f}, Fmscmin: {:.3f}".format(
    #         xmin, Fpcrange[xmin], Ftcsrange[xmin], Fmscrange[xmin]
    #     )
    # )
    xmin = xmin + FMSCMINOFFSET
    if FMSCMINOFFSET != 0:
        print(
            "    OFFSET -> xmin: {:d}, Fpcmin: {:.3f}, Ftcsmin: {:.3f}, Fmscmin: {:.3f}".format(
                xmin, Fpcrange[xmin], Ftcsrange[xmin], Fmscrange[xmin]
            )
        )
    if FIGURE:
        pylab.figure()
        pylab.title("Fpc vs Ftc,s")
        pylab.plot(Fpcrange, Ftcsrange, color='red', marker='o', linestyle='--')
        pylab.plot(Fpcrange[xmin], Ftcsrange[xmin], color='green', marker='x')
        pylab.xlabel('Fpc')
        pylab.ylabel('Ftc,s')
        pylab.grid()
        pylab.legend([f's = {s}, Fps = {round(Fps_fix,5)}', 'minimum'])
        pylab.show()


    N = math.pow(2, NSPACESIZE)  # N value of cipher
    Mchosen = math.pow(2, CHOSENMEMORY)  # chosen memory size tradeoff
    Tchosen = Ftcsrange[xmin] * (N ** 2) / ((Mchosen) ** 2)

    # print("    Memory: {:.3f}, 2^{:.2f}".format(Mchosen, math.log(Mchosen, 2)))
    # print("    Time: {:.3f}, 2^{:.2f}".format(Tchosen, math.log(Tchosen, 2)))

    mm = smallm(Fmscrange[xmin], Fcr_fixedmsc(Fmscrange[xmin]), Ftcsrange[xmin], s, Fps_fix, N, Tchosen)
    tt = smallt(Fps_fix, Fcr_fixedmsc(Fmscrange[xmin]), Ftcsrange[xmin], s, Tchosen)
    ll = smalll(Fps_fix, Fmscrange[xmin], Fcr_fixedmsc(Fmscrange[xmin]), Ftcsrange[xmin], s, Tchosen)
    # print("m: {:.3f}, 2**{:.3f}".format(mm, math.log(mm, 2)))
    # print("t: {:.3f}, 2**{:.3f}".format(tt, math.log(tt, 2)))
    # print("l:{:.3f}, 2**{:.3f}".format(ll, math.log(ll, 2)))
    print('[*] Done!')

    l = ll
    m = mm
    t = tt
    if BLOCKSIZE is False:
        BLOCKSIZE = math.ceil(math.log(N, 2))
    r = int(2 ** (BLOCKSIZE - math.ceil(math.log(t, 2)) - TRUNCATIONBITS))

    EPsize = math.log(r, 2)
    epsilon = math.ceil(math.log(m, 2)) + EPsize - INDEXFILECOMPBITS

    print()
    print("-------- Parameters and estimated performance --------")
    print("m number of chains per table:    2^{:.2f}".format(math.log(m, 2)))
    print("t DP chain length:               2^{:.0f}".format(math.log(t, 2)))
    print("s number of colors:              {:.0f}".format(s))
    print("l number of tables:              2^{:.2f}".format(math.log(l, 2)))
    print("r truncation parameter:          2^{:.0f}".format(math.log(r, 2)))
    print("e endpoint size (bits):          {:.0f}".format(EPsize - INDEXFILECOMPBITS))
    print("------------------------------------------------------")
    print('Available inversion data:        {}'.format(NUMPAIRS))
    # print('Success prob. single:            {:.3f} %'.format(Fps() * 100))
    print('Success rate:                    {:.0f} %'.format(Fps_multi(NUMPAIRS) * 100))
    # print('Fmsc: {:.1f}'.format(Fmsc()))
    print("------------------------------------------------------")
    print('Single table entry size:         {:.0f} bits'.format(epsilon))
    print('Total table size:                {:.2f} TB'.format(m * l * (epsilon) / 1000 ** 4 / 8))
    print('Precomputation cost (Finv^p):    2^{:.2f}'.format(math.log(m * t * l * s, 2)))
    # print('FPGA Enc rate: \t\t\t 2**{:.2f}'.format(math.log(CalcFPGAEncRate(), 2)))
    print('Precomputation FPGA time:        {:.2f} days'.format(FPGATime(m * t * l * s) / 3600 / 24))
    for iterations in (NUMPAIRS, NUMPAIRS * BRUTEFORCE) if BRUTEFORCE > 1 else (NUMPAIRS,):
        print("------------------------------------------------------")
        # print("Number of searches: \t\t {}".format(iterations))
        print('f iter. for online search:       2^{:.2f}'.format(math.log(Tf(iterations), 2)))
        print('f iter. for false alarms:        2^{:.2f}'.format(math.log(Ta(iterations), 2)))
        print('f iter. for truncation alarms:   2^{:.2f}'.format(math.log(Ftrunc(iterations), 2)))
        total_F = Tf(iterations) + Ta(iterations) + Ftrunc(iterations)
        print('Total expected f iter. (Finv^o): 2^{:.2f}'.format(math.log(total_F, 2)))
        print('Total online time FPGA:          {:.2f} s'.format(FPGATime(total_F)))
        print('Total Expected Lookups:          2^{:.2f}'.format(math.log(Lookups(iterations), 2)))
        print('Total Expected Lookups Time:     {:.2f} s'.format(DiskTime(Lookups(iterations)))
        )
        # print('Upper Bound Lookups: \t\t 2**{:.2f}'.format(math.log(LookupsUpperBound(iterations), 2)))
        # print(
        #     'Upper Bound Lookups Time: \t {:.2f}s, {:.2f}m'.format(
        #         DiskTime(LookupsUpperBound(iterations)), DiskTime(LookupsUpperBound(iterations)) / 60
        #     )
        # )
        mintime = max(FPGATime(total_F), DiskTime(Lookups(iterations)))
        maxtime = FPGATime(total_F) + DiskTime(Lookups(iterations))
        mintime_ub = max(FPGATime(total_F), DiskTime(LookupsUpperBound(iterations)))
        maxtime_ub = FPGATime(total_F) + DiskTime(LookupsUpperBound(iterations))
        print('Online Time (T^o):               {:.2f}s <= x <= {:.2f}s'.format(mintime, maxtime))
        # print('Online Wall Clock Time Upper Bound Lookups: {:.2f}s <= x <= {:.2f}s'.format(mintime_ub, maxtime_ub))
        print("------------------------------------------------------")
    print('')
    # compute best s value by minimizing adjusted tradeoff coefficient
    # min_{fmsc}{{log(m_0s)+\varepsilon}^2F_{tc,s}}
    print("[-] Finding optimal s value...")
    m0 = m / s
    srange = range(*SRANGE)  # suggested range of interest 15,101
    Fpcrangemat = []
    Ftcsrangemat = []
    Ftcsminmat = []
    for s in srange:
        print(f'{s}/{srange[-1]}', end='\r')
        Fpcrange = [Fpc_alt(Fps_fix, Fmscel) for Fmscel in Fmscrange]
        Ftcsrange = [
            math.pow((math.log(m0 * s) + EPsize), 2) * Ftc_s_ratio(s, ltratio(Fmscel, Fps_fix), Fmscel)
            for Fmscel in Fmscrange
        ]

        Fpcrangemat.append(Fpcrange)
        Ftcsrangemat.append(Ftcsrange)

        xmin = np.argmin(Ftcsrange)
        Ftcsminmat.append(Ftcsrange[xmin])

    smin = srange[np.argmin(Ftcsminmat)]
    print("[*] optimal s value = ", smin)
    print()
