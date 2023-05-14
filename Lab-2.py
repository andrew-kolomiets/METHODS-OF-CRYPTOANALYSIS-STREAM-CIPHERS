import math
import ctypes

with open('Variants_Gamma_V4.txt', 'r') as gamma_file:
    GAMMA = list(map(int, gamma_file.readline().strip()))

with open('Variants_CombFunction_V4.txt', 'r') as comb_file:
    comb_f = comb_file.readline().replace('} x', '*x').replace('_{', '').replace('}', '').replace('+', '⊕').strip()

print('f = ', comb_f)

def comb_f(x0, x1, x2, x3, x4, x5):
    result = (x0*x3*x5 + x1*x3*x5 + x1*x5 + x2*x3*x5 + x3*x4*x5 + x3*x5 + x3 + x4*x5 + x5 + 1) % 2
    return result

def int_to_bin(num, l):
    return list(map(int,  format(num, f"0{l}b")))

dot = lambda x, y: sum(map(lambda a, b: int(a)*int(b), int_to_bin(x, 6), int_to_bin(y, 6))) % 2

coefs = []
for alpha in range(2**6):
    coef_sum = 0
    for x in range(2**6):
        binary_x = int_to_bin(x, 6)
        comb_f_output = comb_f(*binary_x)
        dot_output = dot(alpha, x)
        sum_output = (comb_f_output + dot_output) % 2
        coef_sum += (-1) ** sum_output
    coefs.append(coef_sum)

print(coefs)

sum_of_squares = sum(c ** 2 for c in coefs)
assert sum_of_squares == 2 ** (2*6)


non_zero_coefs_dict = {}
for i in range(2**6):
    if coefs[i] != 0:
        non_zero_coefs_dict[i] = coefs[i]


print('\n'.join([f'α = {alpha}, W_f(α) = {non_zero_coefs_dict[alpha]}' for alpha in non_zero_coefs_dict]))

weights = [bin(i).count('1') for i in range(2**6)]

sorted_weigths, sorted_coefs = zip(*sorted(zip(weights, coefs)))

i = 0
while i < 2**6:
    if sorted_coefs[i] != 0:
        print('cor(f) =', sorted_weigths[i] - 1)
        break
    i += 1


j = 1
for alpha, coef in non_zero_coefs_dict.items():
    c = 1 if coef < 0 else 0
    g = ''.join([f'x{i} ⊕ ' for i in range(6) if int_to_bin(alpha, 6)[i] == 1])[:-2]
    if c == 1:
        g += f'⊕ {c}'
    prob = (1 + (-1)**c*coef/2**6)/2
    print(f'{j}) g = {g}, Pr{{f(x) = g(x)}} = {prob}')
    j += 1

#в Python даний код не хоче запускатися, використовуємо Jupither і знаходимо polynomials_int
# R.<x> = GF(2^11, 'a')[]
# polynomials = [x^10 + x^3 + 1,
#                 x^10 + x^7 + 1, 
#                 x^10 + x^4 + x^3 + x^1 + 1, 
#                 x^10 + x^8 + x^3 + x^2 + 1,
#                 x^10 + x^8 + x^4 + x^3 + 1,
#                 x^10 + x^8 + x^5 + x^1 + 1]

# polynomials_int_ = [''.join(str(i) for i in poly.list()[:-1]) for poly in polynomials]
# polynomials_int = [int(''.join(str(i) for i in poly.list()[:-1]), 2) for poly in polynomials]

polynomials_int=[576, 516, 864, 706, 610, 786] 

LENGTH = 10

def gen_sequence(poly: int, init_state: int, nbits: int):
    curr_state = init_state
    res_array = (ctypes.c_int * nbits)()

    for i in range(nbits):
        res_array[i] = (curr_state >> (LENGTH - 1)) & 1
        s = 0
        for j in range(LENGTH):
            s ^= ((curr_state & poly) >> j) & 1
        curr_state = (curr_state << 1) | s
    
    return list(res_array)

def hamming_dist(list1, list2):
    dist = 0
    i = 0
    while i < len(list1):
        dist += abs(list1[i] - list2[i])
        i += 1
    return dist

def g1(x3):
    return (x3 + 1) % 2

def g2(x3, x5):
    return (x3 + x5 + 1) % 2

def g3(x1,x4):
    return (x1 + x4) % 2


def g7(x0, x2):
    return (x0 + x2) % 2

T = [int(8 * (0.5 ** (-2)) * math.log((2 ** (10 * i)) / 0.01)) for i in [1, 2]]

print(T)

T1, T2= T

poss_states = {}
for i in range(2**10):
    seq_x3 = gen_sequence(polynomials_int[3], i, T1)
    seq = [g1(seq_x3[j]) for j in range(T1)]
    poss_states[i] = hamming_dist(seq, GAMMA[:T1])

dist = min(poss_states.values())
print('min hamming distance = ', dist/T1)
poss_states = [i for i in poss_states if poss_states[i] == dist]
for s in poss_states:
    print('[X0, X1, X2, X3, X4, X5]: ', f'[-, -, -, -,{s}, -]')



seq_x3 = gen_sequence(polynomials_int[3], 537, T2)
poss_states = {}

for i in range(2**10):
    seq_x5 = gen_sequence(polynomials_int[5], i, T2)
    seq = [g2(seq_x3[j], seq_x5[j]) for j in range(T2)]
    poss_states[i] = hamming_dist(seq, GAMMA[:T2])
dist = min(poss_states.values())
print('min hamming distance = ', dist/T2)
poss_states = [i for i in poss_states if poss_states[i] == dist]
for s in poss_states:
    print('[X0, X1, X2, X3, X4, X5]: ', f'[-, -, -, -, 537, {s}]')

poss_states = {}

for i in range(2**10):
    seq_x1 = gen_sequence(polynomials_int[1], i, T2)
    for j in range(2**10):
        seq_x4 = gen_sequence(polynomials_int[4], j, T2)
        seq = [g3(seq_x1[k], seq_x4[k]) for k in range(T2)]
        poss_states[(i, j)] = hamming_dist(seq, GAMMA[:T2])
dist = min(poss_states.values())
print('min hamming distance = ', dist/T2)
poss_states = [(i, j) for (i, j) in poss_states if poss_states[(i, j)] == dist]
for i, j in poss_states:
    print('[X0, X1, X2, X3, X4, X5]: ', f'[-, {i},-, 537,{j},714]')  
    

seq_x1 = gen_sequence(polynomials_int[1], 764, 100)
seq_x3 = gen_sequence(polynomials_int[3], 537, 100)
seq_x4 = gen_sequence(polynomials_int[4], 996, 100)
seq_x5 = gen_sequence(polynomials_int[5], 714 , 100)

triger = False
for i in range(2**10):
    if triger:
        break
    seq_x0 = gen_sequence(polynomials_int[0], i, 100)
    for j in range(2**10):
        seq_x2 = gen_sequence(polynomials_int[2], j, 100)
        seq = [comb_f(seq_x0[k], seq_x1[k], seq_x2[k], seq_x3[k], seq_x4[k], seq_x5[k]) for k in range(100)]
        if seq == GAMMA[:100]:
            print('[X0, X1, X2, X3, X4, X5]: ', f' [{i}, 764,{j}, 537,996,714]')
            triger = True
            break

gamma_terms = 1000
exponents = [876, 764, 697, 537, 996, 714]

gamma_seq = []
for k in range(gamma_terms):
    terms = []
    for i, exp in enumerate(exponents):
        term = gen_sequence(polynomials_int[i], exp, 1000)[k]
        terms.append(term)
    gamma_seq.append(comb_f(*terms))

if gamma_seq == GAMMA[:gamma_terms]:
    print("True")


