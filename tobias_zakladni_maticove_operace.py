# Základní maticové operace
# Marek Tobiáš, 1. ročník, Bioinformatika
# zimní semestr 2021/22
# Programování NPRG030

from fractions import Fraction

letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
dict, count = {}, 0   # pro ukládání Matrix objektů pod jednotlivými písmeny z letters; počet tvořených matic


class Matrix:   # definuji objekt matice

    def __init__(self, rows, columns, matrix, name = None):
        self.m = rows
        self.n = columns
        self.matrix = matrix    # seznam seznamů jako 2D reprezentace matice
        self.name = name


class Operations:   # vyčlením samostatnou třídu pro operace

    def addition(self, a, b, s):    # sčítání (a=první sčítanec, b=druhý, s=znaménko + nebo -)
        if not self._check(a, b):   # zkontroluji pomocnou funkcí, jestli jsou obě složky matice
            if None in (a, b):  return  # pokud je alespoň jeden ze sčítanců None, vrátím None
            elif isinstance(a, (int, float)) and isinstance(b, (int, float)):   return a+b # pokud jsou oba skaláry vrátíme jejich součet
            else:   return          # v ostatních případech také vrátím None

        if a.m == b.m and a.n == b.n:   # zkontroluji, jestli mají matice vhodné rozměry
            c = [line.copy() for line in a.matrix]  # zkopíruji seznamovou reprezentaci první matice
            for i in range(a.m):        # sečtu každou souřadici
                for j in range(a.n):
                    if s == "+":
                        c[i][j] += b.matrix[i][j]
                    else:
                        c[i][j] -= b.matrix[i][j]

            return Matrix(a.m, a.n, c)  # vrátím objekt Matrix beze jména jako mezivýpočet se spočteným seznamem c
        else:
            print("Error when adding")  # vrátím pokud nemají shodné rozměry


    def multiplication(self, a, b): # násobení
        if not self._check(a, b):
            if None in (a, b):  return
            elif isinstance(a, (int, float)) and isinstance(b, (int, float)):   return a*b
            else:
                alpha, a = self._scalar(a, b)   # pomocnou funkcí zjistím, který z činitelů je skalár a který maticí
                return Matrix(a.m, a.n, [[alpha*element for element in row] for row in a.matrix])   # vrátím matici, která se od původní liší tím, že každý její prvek je vynásoben skalárem alpha

        if a.n == b.m:  # zkontroluji, jestli mají matice vhodné rozměry
            c = []
            for m in range(a.m):
                c.append([])
                for n in range(b.n):
                    sum = 0
                    for j in range(a.n):
                        sum += a.matrix[m][j] * b.matrix[j][n]
                    c[m].append(sum)    # dle definice násobení matic postupně přidávám spočtené sumy do řádků

            return Matrix(a.m, b.n, c)
        else:
            print("Error when multiplying")


    def ref(self, a, i=0, j=0, r=0, expand=1, lu=0): # převod do ref tvaru (a=matice, i,j=začátek na souřadnici A_i,j, r=průběžný rank, expand=sleduji, jestli je matice v rozšířeném tvaru (tj. pokud má hodnotu 2), lu=sleduji, jestli je na matici puštěna funkce lu)
        if not self._check(a):  return None, None
        c = [line.copy() for line in a.matrix]

        while i < a.m and j < a.n:  # dokud nedojdeme na 'konec' matice
            j = self._skip(c, i, j, a.m, a.n // expand) # přeskočíme nulové podsloupečky a najdeme sloupec, ve kterým je pivot
            if j != -1: # pokud jsme v pomocné funkci _skip došli na konec matice, výpočet končí
                for k in range(i+1, a.m):   # s pomocnou funkcí _reduce vynuluji prvky pod pivotem
                    self._reduce(c, i, j, k, a.n, lu)
                i, j, r = i+1, j+1, r+1     # přesunu se o jeden prvek vpravo a jeden dolu a zvýším průběžný rank
            else:
                break

        return Matrix(a.m, a.n, c), r   # vrátím tuple matice a její rank


    def rref(self, a, i=0, j=0, inREF=0):   # převod do rref tvaru (a=matice, i,j=začátek, inREF=pro případ, kdy víme, že je matice už v redukovaném tvaru)
        if not self._check(a):  return
        c = [line.copy() for line in a.matrix]

        while i < a.m and j < a.n:
            j = self._skip(c, i, j, a.m, a.n)
            if j != -1:
                pivot = c[i][j] # najdu pivot a upravím ho na jedničku s pomocí legálních operací
                for x in range(j, a.n):
                    c[i][x] /= pivot
                if not inREF:   # zkontroluji, jestli už náhodou matice není v REF tvaru
                    for k in range(0, a.m): # vynuluji prvky nad i pod pivoty
                        self._reduce(c, i, j, k, a.n)
                else:
                    for k in range(0, i):   # pokud je v REF, ušetřím čas - nuluji pouze nad pivoty
                        self._reduce(c, i, j, k, a.n)
                i, j = i+1, j+1
            else:
                break

        return Matrix(a.m, a.n, c)


    def inversion(self, a): # inverz matice
        if self._check(a) and a.m == a.n:   # zkontroluji, jestli je maticí a jestli má vhodné rozměry
            c = [line.copy() for line in a.matrix]

            expanded, rank = self.ref(Matrix(a.m, 2*a.n, self._expand(c, a.n)), expand=2) # pustím ref na rozšířenou matici tak, abych zjistil, že je invertibilní - musí mít tedy plnou hodnost

            if rank == a.n:
                expanded = self.rref(expanded, inREF=1) # pustím rref na redukovanou rozšířenou matici
                return Matrix(a.m, a.n, [l[a.n:] for l in expanded.matrix]) # vrátím druhou půlku rozšířeného tvaru (inverzní matici)

        print("Not invertible")


    def lu_decomposition(self, a):  # LU dekompozice
        if self._check(a) and a.m == a.n:   # LU dekompozici lze pustit pouze na čvercové matice
            c = [line.copy() for line in a.matrix]

            expanded = self.ref(Matrix(a.m, 2*a.n, self._expand(c, a.n)), expand=2, lu=1)[0]    # pustím ref na rozšířenou matici
            lower, upper = [l[a.n:] for l in expanded.matrix], [l[:a.n] for l in expanded.matrix]
                # rozdělím rozšířenou matici na L a U a vrátím je jako tuple
            return Matrix(a.m, a.n, lower), Matrix(a.m, a.n, upper)
        else:
            print("Error when decomposing!")


    def transposition(self, a): # transpozice
        if not self._check(a):  return a    # definuji transpozici skaláru jako ten samý skalár a pokud je hodnota None, vrátí ji také
        c = []
        for m in range(a.n):    # z řádků sloupce a ze sloupců řádky - přehodím podle diagonály
            c.append([])
            for n in range(a.m):
                c[m].append(a.matrix[n][m])

        return Matrix(a.n, a.m, c)


    # pomocné funkce
    def _skip(self, c, k, l, m, n): # (c=matice, k=nynější řádek , l=nynější sloupec , m=počet řádků, n=počet sloupců)
        i = k
        while c[k][l] == 0: # dokud se prvek rovná 0, přeskakuji
            if k+1 < m: # nejdřív zkontroluji prvky pod počátečním prvkem
                k += 1
            else:
                if l+1 < n: # pokud jsou všechny prvky pod počátečním nulové, přeskočím o sloupec vpravo (pokud to jde)
                    k, l = i, l+1
                else:
                    return -1
        if i != k:  # pokud jsem se přeskakováním dostal na jiný řádek než na kterém jsem začínal, vyměním jej s ním
            c[k], c[i] = c[i], c[k]
        return l    # vrátím sloupec, ve kterém jsem našel pivota

    def _expand(self, c, n):    # rozšířím seznam, který reprezentuje matici, o jednotkovou matici
        for i in range(n):
            for j in range(n):
                if i == j:
                    c[i].append(1)
                else:
                    c[i].append(0)
        return c

    def _reduce(self, c, i, j, k, n, lu=0): # (c=seznam reprezentující matici, i=počáteční řádek, j=počáteční sloupec, k=jaký řádek redukuji, n=šířka matice, lu=jestli je na matici puštěna operace LU dekompozice)
        if k != i:
            alpha = c[k][j] / c[i][j]
            for l in range(j, n):
                if lu and l >= n//2:    # při LU dekompozici horní trojúhelníkové matice
                    if l >= n//2 + j:
                        c[k][l] += alpha*c[i][l]
                else:
                    c[k][l] -= alpha*c[i][l]

    def _check(self, a, b=Matrix(0,0,[])):  # zkontroluji, jestli jsou oba parametry maticemi
        if isinstance(a, Matrix) and isinstance(b, Matrix):
            return True
        else:
            return False

    def _scalar(self, a, b):    # vrátím v tuple - jako první skalár a druhou matici
        if isinstance(a, Matrix):
            return b, a
        else:
            return a, b
    #



class InputLoop:    # samostatná třída pro uživatelský vstup a jeho vyhodnocení

    def __init__(self):
        print("""Type 'matrix' to create one, type 'show' to know what you've created and 'exit' to stop the program.

To apply binary operation, type 'operation(X, Y)'
To apply unary operation, type 'operation(X)'

You can choose from these operations {+, -, *, rank, ref, rref, inverse, lu, transpose}.

Example:
'matrix' ... A = [ 1 2 ]  =>  '+(A, A)'  =>  B = [ 2 4 ]
""")
        self.loop()


    def loop(self): # funkce s while cyklem, který běží dokud exp_eval nevrátí True (program má skončit)
        inProcess = 1
        while inProcess:
            if self.exp_eval(input().replace(" ", "").lower()):
                inProcess = 0


    def create_matrix(self):    # dialog pro tvorbu matice

        def get_dims(): # funkce pro získání dimenzí legálních hodnot
            print("Enter dimensions:")
            m, n = input("m = ").split(), input("n = ").split()

            if len(m) != 1 or len(n) != 1 or any(not dim.isdigit() for dim in [m[0], n[0]]) or int(m[0]) < 1 or int(n[0]) < 1:
                print("Invalid dimension", end="\n\n")
                return get_dims()

            return int(m[0]), int(n[0])

        def get_values():   # získání samotných číselných hodnot každého prvku po řádcích
            c = []
            for i in range(m):
                c.append([])
                for j in range(n):
                    c[i].append("-")
                print("[", *c[i], "]", end=" enter values: ")
                row = input().split()
                length = len(row)
                for x in range(n):
                    if x >= length:
                        c[i][x] = 0
                    else:
                        if row[x].replace("-", "").isdigit():
                            c[i][x] = float(row[x])
                        else:
                            print("Non-numeric value encountered! Enter your matrix again.", end="\n\n")
                            return get_values()

            return c

        m, n = get_dims()

        return Matrix(m, n, get_values())


    def exp_eval(self, exp):    # pro vyhodnocení výrazu

        def command():  # funkce pro zpracování jednotlivých příkazů, které nejsou v rámci operací
            if exp == "matrix": # zahájí dialog pro vytvoření matice
                    self.store([self.create_matrix()])
            elif exp == "show": # ukáže všechny doposud vytvořené matice
                for i in dict:
                    self.visual(dict[i])
            elif exp == "exit" or exp == "end":
                return True
            else:
                print("Unknown command")

        def symbol(interval, sym=",", counter=0):   # najde symbol po korektním uzávorkování
            for i in range(len(interval)):
                if counter < 0:
                    return False
                elif interval[i] == "(":
                    counter += 1
                elif interval[i] == ")":
                    counter -= 1
                elif interval[i] == sym and counter == 0:
                    return i
            if counter != 0:
                return False
            else:
                return True

        def evaluate(i, j, d=0, operator=0):    # (i=začátek intervalu, j=konec intervalu, d=depth, operator="znaménko operace")
            s = i
            while i < j and exp[i] != "(":
                i += 1
            # najdu levou závorku a určím znaménko operace (nebo operaci)
            if i < j:
                operator = exp[s:i]
                i, j = i+1, j-1

            if operator:    # pokud závorku najdu

                if operator == "+" or operator == "-":
                    k = symbol(exp[i:j]) + i    # u binárních operací najdu čárku, která odděluje jednotlivé členy
                    return use.addition(evaluate(i, k, d+1), evaluate(k+1, j, d+1), operator)

                elif operator == "*":
                    k = symbol(exp[i:j]) + i
                    return use.multiplication(evaluate(i, k, d+1), evaluate(k+1, j, d+1))

                elif operator == "rank":
                    return use.ref(evaluate(i, j, d+1))[1]

                elif operator == "ref":
                    return use.ref(evaluate(i, j, d+1))[0]

                elif operator == "rref":
                    return use.rref(evaluate(i, j, d+1))

                elif operator == "inverse":
                    return use.inversion(evaluate(i, j, d+1))

                elif operator == "lu":
                    if d == 0:
                        return use.lu_decomposition(evaluate(i, j, d+1))
                    else:   # nerozkládám pokud se nejedná o poslední operaci ve výrazu
                        return evaluate(i, j, d+1)

                elif operator == "transpose":
                    return use.transposition(evaluate(i, j, d+1))

            else:   # pokud závorku nenajdu -> jde o člen, na kterým se má operace provést
                string = exp[s:j]
                string = string.upper()
                if string in dict:          # pokud je v dict
                    return dict[string]
                elif string.isnumeric():    # pokud je člen skalárem
                    return float(string)
                else:                       # pokud neni ani jedno z toho vrátím None
                    return

        if "(" not in exp:  # pokud není ve výrazu závorka, zkontroluji, o jaký příkaz se jedná
            if command():
                return True
        else:
            if not symbol(exp, " "):    # zkontroluji správné uzávorkování tím, že hledám symbol, který ve výraze není
                return print("Wrong parentheses.")

            result = evaluate(0, len(exp))  # vyhodnotím výraz ve speciální prefixové notaci

            if result == None:
                print("An error occured!")
            elif isinstance(result, Matrix):      # pokud je result maticí, vrátím ho v seznamu (technická pomůcka pro funkci store)
                self.store([result])
            elif isinstance(result, tuple):       # pouze pro případ LU dekompozice
                self.store(result)
            else:
                print(result)                     # pro případ číselného výsledku


    def store(self, result):    # ukládá vzniklé matice tak, že jim přiřadí písmeno z abecedy
        global count
        for mat in result:
            n, count = letters[count % 26], count+1 # po obsazení posledního písmena se začíná znovu prvním
            mat.name, dict[n] = n, mat
            self.visual(dict[n])


    def visual(self, a):    # vizualizuje uložené matice
        maxes = [1]*a.n # najdu maximální délku čísla (převedeného na zlomek, pokud není celé) v každém sloupci
        for col in range(a.n):
            for row in a.matrix:
                length = len(str(Fraction(row[col]).limit_denominator(1000)))
                if length > maxes[col]:
                    maxes[col] = length

        print("")
        for i in range(a.m):
            if i == (a.m-1)//2:
                print(a.name, "= [", end="")
            else:
                print(3*" ", "[", end="")
            for j in range(a.n):    # printnu vycentrované sloupce
                number = str(Fraction(a.matrix[i][j]).limit_denominator(1000))
                spacing, z = (maxes[j] - len(number)) // 2, (maxes[j] - len(number)) % 2
                print((spacing+z)*" ", number, end=spacing*" ") # koriguji proměnnou z
            print(" ]")
        print("")



use = Operations()
user = InputLoop()
