freq = dict()

s = "aaaaaaaabbsdctrrgggggg"

for c in s:
    if c not in freq:
        freq[c] = 1
    else:
        freq[c] += 1

print(freq)


#make a function that detects the frequency of combinations of 2 and 3 letters


def freq2(p):
    freq = {}


    idx = 0

    while idx < len(p) - 1:
        try:
            if p[idx] == p[idx + 1]:
                seq = p[idx] + p[idx]
                if seq not in freq:
                    freq[seq] = 1
                else:
                    freq[seq] += 1

            idx += 1
        except:
            print("finished")
            break

    idx = 0
    while idx < len(p) - 2:
        try:
            # print("here3")
            if p[idx] == p[idx + 1] == p[idx + 2]:
                seq = p[idx] + p[idx] + p[idx]
                if seq not in freq:
                    freq[seq] = 1
                else:
                    freq[seq] += 1
            idx += 1
        except:
            print("finished")
            break
    return freq

print(freq2(s))