import math


def euclidean_distance(X: str, Y: str, k: int) -> float:
    n = len(X)
    m = len(Y)

    if n < k or m < k:
        return float('inf')

    tot_X = n - k + 1
    tot_Y = m - k + 1

    # Estrazione e ordinamento seq X
    w_X = []
    c_X = []
    for i in range(tot_X):
        w = X[i:i + k]
        low = 0
        high = len(w_X)
        pos = low
        while low < high:
            mid = (low + high) // 2
            if w_X[mid] == w:
                pos = mid
                break
            elif w_X[mid] < w:
                low = mid + 1
            else:
                high = mid
        else:
            pos = low

        if pos < len(w_X) and w_X[pos] == w:
            c_X[pos] += 1
        else:
            w_X.insert(pos, w)
            c_X.insert(pos, 1)

    # Estrazione e ordinamento seq Y
    w_Y = []
    c_Y = []
    for j in range(tot_Y):
        w = Y[j:j + k]
        low = 0
        high = len(w_Y)
        pos = low
        while low < high:
            mid = (low + high) // 2
            if w_Y[mid] == w:
                pos = mid
                break
            elif w_Y[mid] < w:
                low = mid + 1
            else:
                high = mid
        else:
            pos = low

        if pos < len(w_Y) and w_Y[pos] == w:
            c_Y[pos] += 1
        else:
            w_Y.insert(pos, w)
            c_Y.insert(pos, 1)

    S = 0.0
    i = 0
    j = 0
    d_X = len(w_X)
    d_Y = len(w_Y)

    while i < d_X and j < d_Y:
        if w_X[i] == w_Y[j]:
            p_X = c_X[i] / tot_X
            p_Y = c_Y[j] / tot_Y
            S += (p_X - p_Y) ** 2
            i += 1
            j += 1
        elif w_X[i] < w_Y[j]:
            p_X = c_X[i] / tot_X
            S += p_X ** 2
            i += 1
        else:
            p_Y = c_Y[j] / tot_Y
            S += p_Y ** 2
            j += 1

    while i < d_X:
        p_X = c_X[i] / tot_X
        S += p_X ** 2
        i += 1

    while j < d_Y:
        p_Y = c_Y[j] / tot_Y
        S += p_Y ** 2
        j += 1

    return math.sqrt(S)


def d2_distance(X: str, Y: str, k: int) -> float:
    n = len(X)
    m = len(Y)

    if n < k or m < k:
        return 1.0

    tot_X = n - k + 1
    tot_Y = m - k + 1

    w_X = [""] * tot_X
    c_X = [0] * tot_X
    d_X = 0

    for i in range(tot_X):
        w = X[i:i + k]
        p = 0
        for j in range(d_X):
            if w_X[j] == w:
                c_X[j] += 1
                p = 1
                break
        if p == 0:
            w_X[d_X] = w
            c_X[d_X] = 1
            d_X += 1

    w_Y = [""] * tot_Y
    c_Y = [0] * tot_Y
    d_Y = 0

    for i in range(tot_Y):
        w = Y[i:i + k]
        p = 0
        for j in range(d_Y):
            if w_Y[j] == w:
                c_Y[j] += 1
                p = 1
                break
        if p == 0:
            w_Y[d_Y] = w
            c_Y[d_Y] = 1
            d_Y += 1

    D2 = 0.0
    sum_sq_X = 0.0
    sum_sq_Y = 0.0

    for i in range(d_X):
        sum_sq_X += c_X[i] ** 2
    for j in range(d_Y):
        sum_sq_Y += c_Y[j] ** 2

    for i in range(d_X):
        w = w_X[i]
        for j in range(d_Y):
            if w_Y[j] == w:
                D2 += c_X[i] * c_Y[j]
                break

    S_X = math.sqrt(sum_sq_X)
    S_Y = math.sqrt(sum_sq_Y)

    if S_X == 0 or S_Y == 0:
        return 1.0

    C = D2 / (S_X * S_Y)
    if C > 1.0:
        C = 1.0

    return 1.0 - C


def shannon_entropy(X: str, Y: str, k: int) -> float:
    n = len(X)
    m = len(Y)

    if n < k or m < k:
        return float('inf')

    tot_X = n - k + 1
    tot_Y = m - k + 1
    U = tot_X + tot_Y

    w_X = [""] * tot_X
    c_X = [0] * tot_X
    d_X = 0
    for i in range(tot_X):
        w = X[i:i + k]
        p = 0
        for j in range(d_X):
            if w_X[j] == w:
                c_X[j] += 1
                p = 1
                break
        if p == 0:
            w_X[d_X] = w
            c_X[d_X] = 1
            d_X += 1

    w_Y = [""] * tot_Y
    c_Y = [0] * tot_Y
    d_Y = 0
    for i in range(tot_Y):
        w = Y[i:i + k]
        p = 0
        for j in range(d_Y):
            if w_Y[j] == w:
                c_Y[j] += 1
                p = 1
                break
        if p == 0:
            w_Y[d_Y] = w
            c_Y[d_Y] = 1
            d_Y += 1

    HX = 0.0
    HY = 0.0
    HU = 0.0

    for i in range(d_X):
        w = w_X[i]
        cnt_X = c_X[i]

        cnt_Y = 0
        for j in range(d_Y):
            if w_Y[j] == w:
                cnt_Y = c_Y[j]
                break

        PX = cnt_X / tot_X
        HX -= PX * math.log2(PX)

        P = (cnt_X + cnt_Y) / U
        HU -= P * math.log2(P)

    for j in range(d_Y):
        w = w_Y[j]
        cnt_Y = c_Y[j]

        PY = cnt_Y / tot_Y
        HY -= PY * math.log2(PY)

        p = 0
        for i in range(d_X):
            if w_X[i] == w:
                p = 1
                break
        if p == 0:
            P = cnt_Y / U
            HU -= P * math.log2(P)

    return HU - 0.5 * (HX + HY)


def average_common_substring(X: str, Y: str) -> float:
    n = len(X)
    m = len(Y)

    if n == 0 or m == 0:
        return float('inf')

    sum_len = 0

    for i in range(n):
        max_len = 0
        for j in range(m):
            h = 0
            while (i + h < n) and (j + h < m) and (X[i + h] == Y[j + h]):
                h += 1
            if h > max_len:
                max_len = h
        sum_len += max_len

    L = sum_len / n
    if L == 0:
        return float('inf')

    return (math.log2(m) / L) - (2 * math.log2(4) / L)


def normalized_compression_distance(X: str, Y: str) -> float:
    def _lz_complexity(S: str) -> int:
        length = len(S)
        M = [""] * length
        dim_M = 0
        i = 0
        cnt = 0

        while i < length:
            h = 1
            while i + h <= length:
                w = S[i:i + h]
                p = 0
                for j in range(dim_M):
                    if M[j] == w:
                        p = 1
                        break
                if p == 0:
                    M[dim_M] = w
                    dim_M += 1
                    break
                h += 1
            cnt += 1
            i += h

        return cnt

    K_X = _lz_complexity(X)
    K_Y = _lz_complexity(Y)
    K_XY = _lz_complexity(X + Y)

    if K_X > K_Y:
        A = K_X
        B = K_Y
    else:
        A = K_Y
        B = K_X

    if A == 0:
        return 1.0

    return (K_XY - B) / A
