import numpy as np
import re


def get_and_conv_user_input(eqIn=None):
    def clean_eq(eqs):
        import re

        ansList = []
        eqs = eqs.replace(" ", "")
        splited = re.split(r'[=><]+', eqs)
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

    # get initial user input
    equations=[]
    print("formatting: all variables on left except for target function where max (min) should be on the left \n"
          "for ex.:\n  min(z)=x1*3+a\n  3*x1-3*a>=100\n  7*a+x2>=35\n  8*a+x2=43\n x1,x2>=0\n"
          "the program will standardize by itself. when done, press enter to run program\n******************")
    while eqIn!="":
        eqIn=input()
        if eqIn=="d":
            equations=[clean_eq(eqIn) for eqIn in ["3*x1-3*a>=100","min(z)=x1*3+a","7*a+x2>=35","8*a+x2=43","x1,x2>=0"]]
            #["max(z)=14*x1+18*x2", "-x1+3*x2<=6", "7*x1+x2<=35"]]
            break
        try: equations.append(clean_eq(eqIn))
        except: print("Inavild Input, try again")

    log_pretabels = True if simp_type=="debug" else False
    if log_pretabels:
        print(*equations, sep="\n")
        print("*--------------------------------*")

    # convert equ to dict to be used later
    flag=False;is_min=False
    for i,a in enumerate(equations):
        if "max" in a["signOfEq"]: equations[i]["signOfEq"]=">";flag=True
        elif "min" in a["signOfEq"]: equations[i]["signOfEq"]="<";flag=True;is_min=True
        equations[i] = standardize_eq(a)
        if flag: equations[i]["signOfEq"]="max"; flag=False

    # place z equ on top
    for i,a in enumerate(equations):
        if a["signOfEq"]=="max":
            equations.insert(0, equations.pop(i))
            break

    if log_pretabels:
        print(*equations, sep="\n")
        print("*--------------------------------*")

    # list of variable
    variables = []
    for eq in equations:
        for k in eq.keys():
            if k=="signOfEq" or k=="sumOfEq": continue
            if not k in variables: variables.append(k)

    # sign constraint
    sign_constraint = []
    for i, eq in enumerate(list(equations)):
        keys=list(eq.keys())
        if len(keys)==3 and eq["sumOfEq"]==0:
            eq = equations.pop(i)

            if "," in keys[0]: vars = re.split(',| , | ,|, ',keys[0])
            else: vars = [keys[0]]

            for v in vars:
                sign_constraint.append(v)

    return equations, sign_constraint, is_min


def create_table(eqs, is_min):
        all_var = set().union(*eqs)
        all_var.remove('signOfEq')
        all_var.remove('sumOfEq')
        all_var = sorted(all_var)

        # turn eqs into a unified table
        initTabel = [[] for i in range(len(eqs))]
        r_count=0
        global simp_type
        for i, e in enumerate(eqs):
            for v in all_var:
                if not v in e: e[v] = 0.0

            if simp_type!="dual" and e["sumOfEq"]<0 or e["signOfEq"]=="=":
                r_count+=1
                for j, eq in enumerate(eqs): eqs[j]["R"+str(r_count)]=0
                eqs[i]["R"+str(r_count)]=1
                all_var.append("R"+str(r_count))
        for i, e in enumerate(eqs):
            for v in all_var:
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

        r_V=[colN.index("R"+str(i+1)) for i in range(r_count)]
        r_count=min([*r_V, 0])-r_count
        for i, eq in enumerate(finalTabel):
            for r in r_V:
                if eq[r]==1:
                    rowN[i]="R"+str(r-r_count-1)
                    if finalTabel[i][-1]<0:
                        for j, a in enumerate(eq): finalTabel[i][j]*=-1
                        finalTabel[i][r] *= -1

        if is_min: rowN[0]="-z"

        return np.array(finalTabel), rowN, colN


def choose_pivot(type_of_simplex, matrix, *args, **kwargs):
    def choose_pivot_simplex(matrix, reverse=-1, *args, **kwargs):
        _, j = max([(a, idx) for (idx, a) in enumerate(matrix[0, :-1]) if reverse*a > 0])
        _, i = min([(abs(r[-1]/r[j]), idx) for (idx, r) in enumerate(matrix[1:]) if r[j]>0 and r[-1]>0])
        return i+1, j

    def choose_pivot_dual(matrix, *args, **kwargs):
        i = min([(a, idx) for (idx, a) in enumerate(matrix[1:,-1]) if a < 0])[1]+1
        _, j = min([(abs(c/matrix[i][idx]), idx) for (idx, c) in enumerate(matrix[0][:-1]) if matrix[i][idx]<0])
        # sort_i = sorted([(a, idx) for (idx, a) in enumerate(matrix[1:,-1]) if a < 0])
        # for val_idx, (a, i) in enumerate(sort_i):
        #     if a != sort_i[0][0]: break
        #     try:
        #         _, j = min([(abs(c/matrix[i][idx]), idx) for (idx, c) in enumerate(matrix[0][:-1]) if c<0])
        #         break
        #     except:
        #         if val_idx==len(sort_i)-1: raise("failed to choose j")
        return i, j
    try:
        return eval("choose_pivot_"+type_of_simplex+"(matrix, *args, **kwargs)")
    except: return None,None


def matrix_reduction(matrix, ii, jj, n_v=False):
    if matrix[ii][jj]==0: print(ii,jj)
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


def print_table(tabel, rowN, colN, ii=None,jj=None, lnn=9):
    print("********************")
    printableTable=[["_", *colN]]
    for i, row in enumerate(tabel):
        printableTable.append([rowN[i], *row])

    row_format = ("{:>"+str(lnn)+"}") * len(printableTable[0])
    print(row_format.format(*printableTable[0]))
    for i, row in enumerate(printableTable[1:]):
        rr=[round(r,3) for r in row[1:]]
        rr.insert(0, row[0])
        if i==ii:
            rr[jj+1]=" "*(lnn-3) +"\033[1m" + str(rr[jj+1]) + "\033[0m"
        print(row_format.format(*rr))


def print_final_var(matrix, colN, rowN):
    print("********************")
    n_v=len(matrix[0])-1
    final_var_list={var:0 for var in colN[:n_v]+[rowN[0]]}
    for i,e in enumerate(matrix[:,-1]):
        final_var_list[rowN[i]]=e
    if "-z" in final_var_list.keys():
        final_var_list["z"]=final_var_list["-z"]
        final_var_list.pop("-z")
    print(final_var_list)


def print_dual(equations, sign_cons):
    print("\n********* Dual Equations *********")

    all_var = set().union(*equations)
    all_var.remove('signOfEq')
    all_var.remove('sumOfEq')
    all_var = sorted(all_var)

    b=[eq["sumOfEq"] for eq in equations[1:]]
    transp_a=[]
    c=[]
    for v in all_var:
        transp_a.append([])
        c.append(equations[0][v])
        for eq in equations[1:]:
            transp_a[-1].append(eq[v])

    dual_eqs=[{}]
    all_var_y=["y"+str(i+1) for i in range(len(transp_a[0]))]

    for i, v in enumerate(all_var_y): dual_eqs[0][v] = b[i]
    dual_eqs[0]["signOfEq"]="min y0"
    dual_eqs[0]["sumOfEq"]=0
    for idx, eq in enumerate(transp_a):
        new_eq = {}
        for i, v in enumerate(all_var_y): new_eq[v] = eq[i]
        new_eq["signOfEq"]=">="
        new_eq["sumOfEq"]=c[idx]
        dual_eqs.append(new_eq)


    print("min y0=", end="")
    for k, v in zip(list(dual_eqs[0])[:-2],list(dual_eqs[0].values())[:-2]):
        if v==1: v="+"
        elif v==-1: v="-"
        print(str(v)+"*"+k, end="")
    print("")
    for i, eq in enumerate(dual_eqs[1:]):
        for k, v in zip(list(eq)[:-2],list(eq.values())[:-2]):
            if v==1: v="+"
            elif v==-1: v="-"
            elif v>0: v="+"+str(v)+"*"
            elif v<0: v=str(v)+"*"
            elif v==0: continue
            print(str(v)+k, end="")
        if all_var[i] in sign_cons: print(">="+str(eq["sumOfEq"]))
        else: print("="+str(eq["sumOfEq"]))

    for i, eq in enumerate(equations):
        if eq["signOfEq"]=="=": print("y"+str(i)+" is free in sign")


def phase_one_twoPhase(table, colNames, rowNames):
    from re import match
    w_0 = [-1 if match("R[1-]", v) else 0 for v in colNames]
    table = np.array([np.array(w_0), *table])
    rowNames = ["w0"] + rowNames

    ii, jj = choose_pivot(simp_type, table)
    print_table(table, rowNames, colNames, ii, jj)
    table = matrix_reduction(table, ii, jj)
    rowNames[ii] = colNames[jj]
    loop_flag = True
    while loop_flag:
        ii, jj = choose_pivot(simp_type, table, reverse=1)
        print_table(table, rowNames, colNames, ii, jj)
        if ii == None: break
        table = matrix_reduction(table, ii, jj)
        rowNames[ii] = colNames[jj]

    table = table[1:]
    r1 = [i for i, v in enumerate(colNames) if match("R[1-]", v)]
    r2 = min(r1); r1 = min(r1)
    for i in range(r1, r2 + 1): table = np.delete(table, r1, 1)
    rowNames.pop(0)
    colNames = [v for v in colNames if not match("R[1-]", v)]
    return table, colNames, rowNames


if __name__=="__main__":
    global simp_type
    simp_type = input("what type of simplex do you want? \n (example: dual, simplex, twoPhase)\n")
    simp_type = {"si": "simplex", "du": "dual", "tw": "twoPhase", "de": "debug"}[simp_type[:2]]
    print("chosen type: ", simp_type)

    equations, sign_constraint, is_min = get_and_conv_user_input()
    table, rowNames, colNames = create_table(equations, is_min)

    if simp_type=="dual": print_dual(equations, sign_constraint)

    if simp_type=="twoPhase":
        simp_type="simplex"
        print("****Phase 1*******")
        table, colNames, rowNames = phase_one_twoPhase(table, colNames, rowNames)
        print("****Phase 2*******")

    loop_flag = True
    while loop_flag:
        ii, jj = choose_pivot(simp_type, table)
        print_table(table, rowNames, colNames, ii, jj)
        if ii == None: break
        table = matrix_reduction(table, ii, jj)
        rowNames[ii] = colNames[jj]

    print_final_var(table, colNames, rowNames)
