def edit_distance(X: str, Y: str) -> list:
    n = len(X)
    m = len(Y)

    M = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):  # inizializza colonna 0
        M[i][0] = i

    for j in range(m + 1):  # inizializza riga 0
        M[0][j] = j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if X[i - 1] == Y[j - 1]:
                p = 0
            else:
                p = 1

            M[i][j] = min(
                M[i - 1][j] + 1,
                M[i - 1][j - 1] + p,
                M[i][j - 1] + 1
            )

    return M


def needleman_wunsch(X: str, Y: str, match=2, mismatch=-1, gap=-1) -> list:
    n = len(X)
    m = len(Y)
    M = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        M[i][0] = i * gap
    for j in range(m + 1):
        M[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if X[i - 1] == Y[j - 1]:
                p = match
            else:
                p = mismatch

            M[i][j] = max(
                M[i - 1][j] + gap,
                M[i - 1][j - 1] + p,
                M[i][j - 1] + gap
            )
    return M


def smith_waterman(X: str, Y: str, match=2, mismatch=-1, gap=-1) -> list:
    n = len(X)
    m = len(Y)
    M = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if X[i - 1] == Y[j - 1]:
                p = match
            else:
                p = mismatch

            M[i][j] = max(
                0,
                M[i - 1][j] + gap,
                M[i - 1][j - 1] + p,
                M[i][j - 1] + gap
            )
    return M


def allinea(X: str, Y: str, M: list, ALX: list, ALY: list) -> int:
    n = len(X)
    m = len(Y)

    k = n + m - 1
    h = k
    i = n
    j = m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and (
                (M[i][j] == M[i - 1][j - 1] and X[i - 1] == Y[j - 1]) or
                (M[i][j] == M[i - 1][j - 1] + 1 and X[i - 1] != Y[j - 1])
        ):
            ALX[h] = X[i - 1]
            ALY[h] = Y[j - 1]
            i = i - 1
            j = j - 1

        else:
            if j > 0 and M[i][j] == M[i][j - 1] + 1:
                ALX[h] = "-"
                ALY[h] = Y[j - 1]
                j = j - 1

            elif i > 0:
                ALX[h] = X[i - 1]
                ALY[h] = "-"
                i = i - 1
        h = h - 1
    return h


def allinea_needleman(X: str, Y: str, M: list, ALX: list, ALY: list, match=2, mismatch=-1, gap=-1) -> int:
    n = len(X)
    m = len(Y)
    k = n + m - 1
    h = k
    i = n
    j = m

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            p = match if X[i - 1] == Y[j - 1] else mismatch
            if M[i][j] == M[i - 1][j - 1] + p:
                ALX[h] = X[i - 1]
                ALY[h] = Y[j - 1]
                i -= 1
                j -= 1
                h -= 1
                continue

        if j > 0 and M[i][j] == M[i][j - 1] + gap:
            ALX[h] = "-"
            ALY[h] = Y[j - 1]
            j -= 1
        elif i > 0:
            ALX[h] = X[i - 1]
            ALY[h] = "-"
            i -= 1
        h -= 1
    return h


def allinea_smith(X: str, Y: str, M: list, ALX: list, ALY: list, match=2, mismatch=-1, gap=-1) -> tuple:
    n = len(X)
    m = len(Y)
    k = n + m - 1
    h = k

    max_val = -1
    start_i, start_j = 0, 0
    for r in range(n + 1):
        for c in range(m + 1):
            if M[r][c] > max_val:
                max_val = M[r][c]
                start_i, start_j = r, c

    i, j = start_i, start_j

    while (i > 0 or j > 0) and M[i][j] != 0:
        if i > 0 and j > 0:
            p = match if X[i - 1] == Y[j - 1] else mismatch
            if M[i][j] == M[i - 1][j - 1] + p:
                ALX[h] = X[i - 1]
                ALY[h] = Y[j - 1]
                i -= 1
                j -= 1
                h -= 1
                continue

        if j > 0 and M[i][j] == M[i][j - 1] + gap:
            ALX[h] = "-"
            ALY[h] = Y[j - 1]
            j -= 1
        elif i > 0:
            ALX[h] = X[i - 1]
            ALY[h] = "-"
            i -= 1
        h -= 1

    return h, i, j
