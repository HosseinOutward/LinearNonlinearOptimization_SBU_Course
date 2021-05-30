def get_and_conv_user_input(eqIn=None, log_pretabels=False):
    def clean_eq(eqs):
        import parser
        import re

        ansList = []
        eqs = eqs.replace(" ", "")
        splited = re.split(r'[=><]+', eqs)
        varNames = parser.expr(splited[0]).compile().co_names
        finalDict = {}

        ansList.append(eqs[len(splited[0]):-len(splited[1])])
        if "min" in splited[0] or "max" in splited[0]:
            ansList[0]=splited[0]
            ansList.append(0)
            splited[0]=splited[-1]
        else: ansList.append(float(splited[-1]))

        splited = splited[0].split("+")
        for s in splited:
            if "-" in s:
                for a in s[s.index("-")+1:]:
                    if not a.isnumeric(): break
                if s.index("-")==0 and a!="*":
                    ansList.insert(0, s.replace("-", "-1*"))
                    continue
                if a!="*":
                    s = s.replace("-", "+-1*").split("+")
                else:
                    s=s.replace("-","+-").split("+")
                ansList.insert(0, s.pop())
                ansList.insert(0, s.pop())

            else: ansList.insert(0, s)

        for s in ansList[:-2]:
            s = s.replace("(", "").replace(")", "")
            if "*" in s:
                s = s.split("*")
                try:
                    finalDict[s[1]] = float(s[0])
                except:
                    finalDict[s[0]] = float(s[1])
            else:
                finalDict[s]=1

        finalDict["sumOfEq"] = ansList[-1]
        finalDict["signOfEq"] = ansList[-2]

        return finalDict

    def standardize_eq(eqs):
        if ">" in eqs["signOfEq"]:
            eqs["signOfEq"] = eqs["signOfEq"].replace(">", "<")
            for a in list(eqs.keys())[:-1]:
                eqs[a] = -eqs[a]
        return eqs

    def create_table(eqs):
        all_var = set().union(*eqs)
        all_var.remove('signOfEq')
        all_var.remove('sumOfEq')
        all_var = sorted(all_var)

        # turn eqs into a unified tabel
        initTabel = [[] for i in range(len(eqs))]
        for i, e in enumerate(eqs):
            for v in all_var:
                if not v in e: e[v] = 0.0
                initTabel[i].append(e[v])
            initTabel[i].append(e["sumOfEq"])

        # create a simplex table
        z = initTabel.pop(0)
        rhs = z.pop()
        for j in range(0, len(eqs)-1): z.append(0)
        z.append(rhs)

        finalTabel = []
        for i in range(len(eqs) - 1):
            # reorder table to be solvable
            for j, e in enumerate(initTabel):
                if e[i] != 0:
                    finalTabel.append(initTabel.pop(j))
                    flag = True
                    break
                elif e[:i].count(0) != len(e):
                    for k, fe in enumerate(finalTabel):
                        if fe[i]!=0 and e[k]!=0:
                            flag = True
                            fe=finalTabel.pop(k)
                            finalTabel.insert(k, initTabel.pop(j))
                            finalTabel.append(fe)
                            break
            if flag is False and z[i] != 0:
                raise ValueError("is not linearly independent (" + all_var[i] + " has no constraint)")

        for i in range(len(finalTabel)):
            # add s0 s1,... to table
            rhs = finalTabel[i].pop()
            for j in range(0, i): finalTabel[i].append(0.0)
            finalTabel[i].append(1.0)
            for j in range(i + 1, len(eqs) - 1): finalTabel[i].append(0.0)
            finalTabel[i].append(rhs)
        finalTabel.insert(0, z)

        colN = list(all_var)
        for i in range(len(eqs) - 1): colN.append("s" + str(i + 1))
        colN.append("RHS")
        rowN=["z"]
        for i in range(len(eqs) - 1): rowN.append("s" + str(i + 1))

        return finalTabel, rowN, colN

    # get initial user input
    equations=[]
    print("formatting: all variables on left except for target function where max (min) should be on the left \n"
          "for ex.:\n  min(z)=x1*3+a \n  3*x1-3*a>=100\n  7*a+x2>=35\n "
          "the program will standardize by itself. when done, press enter to run program\n******************")
    while eqIn!="":
        if eqIn=="d":
            equations=[clean_eq(eqIn) for eqIn in ["min(z)=x1*3+a","3*x1-3*a>=100","7*a+x2>=35"]]
            #["max(z)=14*x1+18*x2", "-x1+3*x2<=6", "7*x1+x2<=35"]]
            break
        eqIn=input()
        try: equations.append(clean_eq(eqIn))
        except: print("Inavild Input, try again")

    if log_pretabels:
        print(*equations, sep="\n")
        print("*--------------------------------*")

    # convert equ to dict to be used laterx
    for i,a in enumerate(equations):
        flag=False
        if "max" in a["signOfEq"]: equations[i]["signOfEq"]=">";flag=True
        elif "min" in a["signOfEq"]: equations[i]["signOfEq"]="<";flag=True
        equations[i] = standardize_eq(a)
        if flag: equations[i]["signOfEq"]="max"

    # place z equ on top
    for i,a in enumerate(equations):
        if a["signOfEq"]=="max":
            equations.insert(0, equations.pop(i))
            break

    if log_pretabels:
        print(*equations, sep="\n")
        print("*--------------------------------*")

    # convort dict to final table
    equations, rowNa, colNa = create_table(equations)

    if len(equations)!=len(rowNa): raise ValueError("insufficient number of equ. for the number of variables")

    return equations, rowNa, colNa


def print_tabel(tabel, rowN, colN,lnn=9):
    print("********************")
    printableTable=[["_", *colN]]
    for i, row in enumerate(tabel):
        printableTable.append([rowN[i], *row])

    row_format = ("{:>"+str(lnn)+"}") * len(printableTable[0])
    print(row_format.format(*printableTable[0]))
    for row in printableTable[1:]:
        rr=[round(r,3) for r in row[1:]]
        rr.insert(0, row[0])
        print(row_format.format(*rr))


def choose_pivot_simplex(matrix, prev_selected):
    p_s_i = [i for i, j in prev_selected]
    _, j = min([(a, idx) for (idx, a) in enumerate(matrix[0][:-1])])
    _, i = min([(e[-1]/e[j], idx) for (idx, e) in enumerate(matrix) if e[j]!=0 and idx!=0 and not idx in p_s_i])
    return i, j


def choose_pivot_dual(matrix, prev_selected, n_v):
    p_s_i = [i for i, j in prev_selected]
    _, i = min([(e[-1],idx) for (idx, e) in enumerate(matrix) if idx!=0 and not idx in p_s_i])
    min_j=float("inf")
    for idx, a in enumerate(matrix[0][:n_v]):
        if a == 0 or matrix[i][idx] == 0 or abs(a / matrix[i][idx]) > min_j: continue
        min_j = abs(a / matrix[i][idx])
        j = idx
    return i, j


def add_u_and_pivot(matrix, prev_selected, rowN, colN, n_v):
    decimal_list=[(round(e[-1]-e[-1]//1,13), idx) for (idx, e) in enumerate(matrix[:n_v+1]) if idx!=0]
    # if
    _, chosen_i = max(decimal_list)

    rowN.append("u"+str(len(prev_selected)+1))
    matrix.append([a//1-a for a in matrix[chosen_i]])

    colN.insert(-1,"u"+str(len(prev_selected)+1))
    for row in matrix:
        rhs=row.pop()
        row.append(0)
        row.append(rhs)
    matrix[-1][-2]=1

    min_j=float("inf")
    n_v+=len(prev_selected)
    for idx, a in enumerate(matrix[0][n_v+1:n_v*2+1], start=n_v+1):
        if a==0 or matrix[-1][idx]==0 or abs(a / matrix[-1][idx])>min_j: continue
        min_j=abs(a / matrix[-1][idx])
        j=idx
    i = len(matrix)-1

    return matrix, rowN, colN, chosen_i, i, j


def matrix_reduction(matrix, ii, jj, n_v=False):
    for i, row in enumerate(matrix):
        if i == ii: continue
        count = matrix[i][jj]/matrix[ii][jj]
        if n_v and 0<i<=n_v and row[-1].is_integer(): continue
        for j, a in enumerate(row):
            matrix[i][j] -= count*matrix[ii][j]

    count=1/matrix[ii][jj]
    for j, a in enumerate(row):
        matrix[ii][j] *= count

    return matrix


def print_final_var(matrix, colN, n_e, n_v):
    print("********************")
    final_var_list={var:0 for var in colN[:n_v]}
    for i,e in enumerate(matrix[1:], start=1):
        for j, a in enumerate(e[:n_v]):
            if matrix[0][j]!=0: continue
            if a==1: final_var_list[colN[j]]=round(matrix[i][-1], 12)
    print(final_var_list)


if __name__=="__main__":
    table, rowNames, colNames = get_and_conv_user_input("d", True)
    number_of_var = (len(colNames) - 1) // 2
    number_of_eq = len([True for a in rowNames if "s" in a])+1
    print_tabel(table, rowNames, colNames)

    loop_flag = True
    while loop_flag:
        loop_flag=False
        pivot = []
        while len([True for a in table[1:] if round(a[-1],13) < 0]) != 0:
            loop_flag=True
            pivot.append(choose_pivot_dual(table, pivot, number_of_var))
            table = matrix_reduction(table, *pivot[-1])
            print_tabel(table, rowNames, colNames)

        pivot = []
        while len([True for a in table[0][:-1] if round(a,13) < 0]) != 0:
            loop_flag=True
            pivot.append(choose_pivot_simplex(table, pivot))
            table = matrix_reduction(table, *pivot[-1])
            print_tabel(table, rowNames, colNames)

        pivot = []
        while len([True for a in table[1:number_of_var+1] if not round(a[-1],13).is_integer()]) != 0:
            loop_flag=True
            table, rowNames, colNames, chosen_i, i, j = \
                add_u_and_pivot(table, pivot, rowNames, colNames, number_of_var)
            print_tabel(table, rowNames, colNames)
            pivot.append((i,j))
            table = matrix_reduction(table, *pivot[-1], number_of_var)
            print_tabel(table, rowNames, colNames)

    print_final_var(table, colNames, number_of_eq, number_of_var+1)