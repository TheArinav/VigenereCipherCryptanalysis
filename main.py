import math
from functools import reduce
from idlelib.help import copy_strip
from importlib.metadata import distribution


def kasiski_examine(ct):
    SEGMENT_SIZE = 3

    # Finds repeating segments of length SEGMENT_SIZE
    def find_reps(_divided):
        _reps = {}
        for i in range(len(_divided)):
            cur = _divided[i]
            for j in range(i + 1, len(_divided)):
                if cur == _divided[j]:
                    if cur in _reps:
                        _reps[cur].append(j * SEGMENT_SIZE)
                    else:
                        _reps[cur] = [i * SEGMENT_SIZE, j * SEGMENT_SIZE]
        return _reps

    # Finds the greatest common divisor of a list of integers
    def gcd_of_vector(vec):
        return reduce(math.gcd, vec)

    # Finds all divisors of a given number
    def find_divisors(n):
        divisors = set()
        for i in range(2, int(math.sqrt(n)) + 1):
            if n % i == 0:
                divisors.add(i)
                divisors.add(n // i)
        if n > 1:
            divisors.add(n)
        return sorted(divisors)

    # Divide the ciphertext into substrings of length SEGMENT_SIZE
    divided = [ct[i:i + SEGMENT_SIZE] for i in range(0, len(ct) - SEGMENT_SIZE + 1)]
    reps = find_reps(divided)
    phrases_to_distance = {}
    klens_by_phrases = []

    # Calculate distances between repeated segments
    for key, positions in reps.items():
        if len(positions) > 1:
            phrases_to_distance[key] = [positions[i] - positions[i - 1] for i in range(1, len(positions))]

    # For each group of distances, calculate the GCD and find divisors
    for distances in phrases_to_distance.values():
        if distances:
            gcd_value = gcd_of_vector(distances)
            if gcd_value > 1:
                klens_by_phrases.append(find_divisors(gcd_value))

    # Find the union of all groups of divisors to consider multiple possible key lengths
    if klens_by_phrases:
        key_lengths = set(klens_by_phrases[0])
        for i in range(1, len(klens_by_phrases)):
            key_lengths.update(klens_by_phrases[i])

        return sorted(key_lengths)
    return []

def freidman_test(ct):
    freq = {}
    N = len(ct)
    IC = 0

    # Count the frequency of each letter in the ciphertext
    for char in ct:
        if char in freq:
            freq[char] += 1
        else:
            freq[char] = 1

    # Calculate the Index of Coincidence (IC)
    for ni in freq.values():
        IC += ni * (ni - 1)

    IC /= (N * (N - 1))  # Normalize the IC by dividing by N(N-1)

    # Estimate key length using Friedman's formula
    key_length_estimate = (0.027 * N) / ((0.065 - IC) + N * (IC - 0.038))

    return int(key_length_estimate)

def sort_by_distance(arr, X):
    return sorted(arr, key=lambda n: abs(n - X))

def key_lengths(ct):
    return sort_by_distance(kasiski_examine(ct),freidman_test(ct))

def is_english(message, word_percentage=20, letter_percentage=85):
    UPPER_LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    LETTERS_AND_SPACE = UPPER_LETTERS + UPPER_LETTERS.lower() + ' \t\n'

    def load_dictionary():
        dictionary_file = open('dictionary.txt')
        english_words = {}
        for word in dictionary_file.read().split('\n'):
            english_words[word] = None
        dictionary_file.close()
        return english_words

    ENGLISH_WORDS = load_dictionary()

    def get_english_count(message):
        message = message.upper()
        message = remove_non_letters(message)
        possible_words = message.split()

        if not possible_words:
            return 0.0  # No words at all, so return 0.0.

        matches = 0
        for word in possible_words:
            if word in ENGLISH_WORDS:
                matches += 1
        return float(matches) / len(possible_words)

    def remove_non_letters(message):
        letters_only = []
        for symbol in message:
            if symbol in LETTERS_AND_SPACE:
                letters_only.append(symbol)
        return ''.join(letters_only)

    words_match = get_english_count(message) * 100 >= word_percentage
    num_letters = len(remove_non_letters(message))
    message_letters_percentage = float(num_letters) / len(message) * 100
    letters_match = message_letters_percentage >= letter_percentage
    return words_match and letters_match

def separate_ceaser(ct: str, k_len: int):
    strips = []
    for i in range(k_len):
        strips.append([])
    for i in range(len(ct)):
        strips[i % k_len].append(ct[i])
    return strips

def get_normal_char_distribution():
    return {
        'E' : 12.7, 'T' : 9.1,  'A' : 8.2, 'O' : 7.5,
        'I' : 7.0,  'N' : 6.7,  'S' : 6.3, 'R' : 6.0, 'H' : 6.1,
        'D' : 4.3,  'L' : 4.0,  'C' : 2.8, 'U' : 2.8, 'M' : 2.4,
        'W' : 2.4, 'F' : 2.2,  'G' : 2.0, 'Y' : 2.0, 'P' : 1.9,
        'B' : 1.5,  'V' : 0.98,  'K' : 0.77, 'X' : 0.15, 'J' : 0.15,
        'Q' : 0.095,  'Z' : 0.074
    }

def get_strip_char_distribution(strip):
    strip = "".join(strip)
    res_freq = {}
    for c in strip:
        if c in res_freq.keys():
            res_freq[c] += 1
        else:
            res_freq[c] = 0
    for k in res_freq:
        res_freq[k] = (res_freq[k] / len(strip)) * 100

    return res_freq

def closest_match(ct_char, ct_char_freq):
    freq = get_normal_char_distribution()
    freq_sorted = [c for c in freq.keys()]
    freq_sorted = sorted(freq_sorted, key=lambda n: (abs(ct_char_freq - freq[n])))
    return freq_sorted

def normalize_frequencies(observed_freqs):
    reg_freq = get_normal_char_distribution()
    target_min = min(reg_freq.values())
    target_max = max(reg_freq.values())
    observed_min = min(observed_freqs.values())
    observed_max = max(observed_freqs.values())
    normalized_freqs = {}
    for (key, value) in dict(observed_freqs).items():
        x_normalized = ((value - observed_min) / (observed_max - observed_min)) * (target_max - target_min) + target_min
        normalized_freqs[key]=x_normalized

    return normalized_freqs

def gen_candidates(strip_freq):
    chars = [c for c in strip_freq.keys()]
    candidates = {}
    normal_sf = normalize_frequencies(strip_freq)
    for c in chars:
        candidates[c] = closest_match(c, normal_sf[c])
    return candidates

def guess_key_char(strip_freq):
    cand = gen_candidates(strip_freq)

# Example usage
ciphertext = '''

J FTNSBCSG BG P DVXMR POR PT O NPICH PJERXOU CBHJSOAJGI TDTORXOU PMZ BZ HXNS DCGTSJXOU POR IFGIJBV UVT XCGMR PSCJOR BFADWWCH DXFQTT OAUSGJBV UVT GZDX CU UVXOUH BBS ECRVATOHXOU LBMH UVT XCGMR GFGEPBSFR IP AT OCL BG PO OSVZI BBS B DGPTTTGXPBPM BPUIGBZXTH XWS PQDGPORISS MOCHIPHS XO HWF GPNS LBM CPH USCB BB PDOSFAXD DDJBI PT KJSL CII BG P DIGJCJT QWJZS THXMZ QVWAEWCH ZXUHAF AJE RPNG XO QGFSZT OCE QWBGXOU PGHTS TGPUH TC IIWH CCDL WH BB DER IIWCH WI JG P OOIVFPMWHUG LBZZ UVGPIVI HWF ZPOUJBUTNOZJBV MOCEGRBDT PT IIS TOUAJGW MOCHIPHS POR UPZAPKXOU XO HWF BPUIGBZXTHH UFPEWIJCC JH RPAQJBTT CQTSGWOIJCC FLEFFXNSCUOIJCC TDTDIABHXPB POR SPQJNSCUOIJCCBQIJJXUWTT KT ECCU BDSAPMZN BGHPQXBHT XWII ZPOUJBUT

'''.replace('\n', '').replace(' ', '')

print(key_lengths(ciphertext))
strips = separate_ceaser(ciphertext,key_lengths(ciphertext)[1])
print("\n\n".join([str(gen_candidates(get_strip_char_distribution(c_strip))).replace('],','],\n') for c_strip in strips]))
