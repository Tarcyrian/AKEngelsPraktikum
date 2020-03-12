import os
import sys
import subprocess
import copy
import numpy as np
import math as m

'''Klasse Atom_xyz:
    Beinhaltet die eingelesenen Symbole
    und Koordinaten in einer Liste'''
class atom_xyz(object):
    def __init__(self, symbol, coords):
        self.symbol = symbol
        self.coords = coords

'''Klasse charge_xyz:
    Beinhaltet die Koordinaten und die
    dazugehörgen Ladungen'''
class charge_xyz(object):
    def __init__(self, coords, charge):
        self.coords = coords
        self.charge = charge

'''Funktion readXYZ:
    Liest die angegebene Datei ein und speichert
    sowohl die Symbole als auch die Koordinaten seperat ab
    und fügt sie anschließend an eine Liste ret an.''' 
def readXYZ(name):
    ret = []
    f = open(name, "r")
    for line in f:
        tmp = line.split()
        if len(tmp) == 4:
            sym = tmp[0]
            coords = list(map(float, tmp[1:]))
            ret.append(atom_xyz(sym, coords))
    f.close()

    return ret

'''Funktion readcharges:
    List die angegebene Datei ein und speichert
    sowohl die Koordinaten als auch die Ladungen seperat ab
    und fügt sie anschließend an eine Liste charges an.'''
def readcharges(file):
    charges = []
    f = open(file, "r")
    for line in f:
        tmp = line.split()
        if len(tmp) == 4:
            coord = list(map(float, tmp[0:3]))
            charge = float(tmp[3])
            charges.append(charge_xyz(coord, charge))
    f.close()

    return charges

'''Funktion getcenterofGeo:
    Berechnet Center of Geometry und gibt die
    veränderten Koordinaten zurück.'''
def getcenterofGeo(ret):
    centerOfGeo = [0.0, 0.0, 0.0]

    for i in range(len(ret)):
        j = 0
        while j < 3:
            centerOfGeo[j] += ret[i].coords[j]
            j += 1

    k = 0       
    while k < 3:
            centerOfGeo[k] /= len(ret)
            k += 1     

    negativeCenterOfGeo = copy.deepcopy(centerOfGeo)

    l = 0       
    while l < 3:
            negativeCenterOfGeo[l] *= -1.0
            l += 1 

    for i in range(len(ret)):
        m = 0
        while m < 3:
            ret[i].coords[m] += negativeCenterOfGeo[m]
            m += 1

    return ret

'''Funktion duplicateDIP:
    Dupliziert ret und verschiebt x oder y um einen konstanten Wert.
    Hängt die veränderten Koordinaten an die ursprüngliche Datei an 
    und gibt sie zurück.'''
def duplicateDIP(ret, verschiebung, x, y):
    ret2 = copy.deepcopy(ret)

    for i in range(len(ret)):

        if verschiebung == "xp":
            ret2[i].coords[0] += x
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "xn":
            ret2[i].coords[0] -= x
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "yp":
            ret2[i].coords[1] += y
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "yn":
            ret2[i].coords[1] -= y
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

    return ret

'''Funktion duplicatePDIR:
    Dupliziert ret und verschiebt x oder y um einen konstanten Wert.
    Hängt die veränderten Koordinaten an die ursprüngliche Datei an 
    und gibt sie zurück.'''
def duplicatePDIR(ret, x, difb, difc, verschiebung):
    ret2 = copy.deepcopy(ret)

    for i in range(len(ret)):
        if verschiebung == "ap":
            ret2[i].coords[0] += x
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "an":
            ret2[i].coords[0] -= x
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "bp":
            ret2[i].coords[0] += difb[0]
            ret2[i].coords[1] += difb[1]
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "bn":
            ret2[i].coords[0] -= difb[0]
            ret2[i].coords[1] -= difb[1]
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "cp":
            ret2[i].coords[0] += difc[0]
            ret2[i].coords[1] += difc[1]
            ret2[i].coords[2] += difc[2]
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

        elif verschiebung == "cn":
            ret2[i].coords[0] -= difc[0]
            ret2[i].coords[1] -= difc[1]
            ret2[i].coords[2] -= difc[2]
            ret.append(atom_xyz(ret2[i].symbol, ret2[i].coords))

    return ret


'''Bildet die Ladungen um DIP1 und speichert sie in einer Liste.'''  
def getchargesDIP1(charge, a, b):
    newcharge =[]
    charge1 = copy.deepcopy(charge)
    charge2 = copy.deepcopy(charge)

    for i in range(len(charge)):
        charge1[i].coords[0] += a
        newcharge.append(charge_xyz(charge1[i].coords, charge1[i].charge))

    for i in range(len(charge)):   
        charge2[i].coords[0] -= a
        newcharge.append(charge_xyz(charge2[i].coords, charge2[i].charge))
       
    return newcharge

'''Bildet die Ladungen um DIP2 und speichert sie in einer Liste.'''
def getchargesDIP2(charge, a, b):
    #charge1 = copy.deepcopy(charge)
    charge2 = copy.deepcopy(charge)
    charge3 = copy.deepcopy(charge)
    #charge4 = copy.deepcopy(charge)
    #charge5 = copy.deepcopy(charge)
    charge6 = copy.deepcopy(charge)

    '''for i in range(len(charge)):
        charge1[i].coords[0] += a
        charge.append(charge_xyz(charge1[i].coords, charge1[i].charge))'''
        
    for i in range(len(charge2)):
        charge2[i].coords[0] -= a
        charge.append(charge_xyz(charge2[i].coords, charge2[i].charge))
        
    for i in range(len(charge3)):
        charge3[i].coords[1] += b
        charge.append(charge_xyz(charge3[i].coords, charge3[i].charge))

    '''for i in range(len(charge4)):
                    charge4[i].coords[1] -= b
                    charge.append(charge_xyz(charge4[i].coords, charge4[i].charge))
                    
                for i in range(len(charge5)):
                    charge5[i].coords[1] -= b
                    charge5[i].coords[0] += a
                    charge.append(charge_xyz(charge5[i].coords, charge5[i].charge))'''

    for i in range(len(charge6)):
        charge6[i].coords[1] += b
        charge6[i].coords[0] -= a
        charge.append(charge_xyz(charge6[i].coords, charge6[i].charge))

    return charge

'''Bildet die Ladungen um PDIR und speichert sie in einer Liste.
Die Ladungen werden anhand der Rotationsaxe ausgewählt.'''
def getchargesPDIR(charge, difb, a, rotAxis):
    newcharge = []
    
    difbplusa = difb[0] + a
    difbminusa = difb[0] - a
    minusdifbminusa = - difb[0] - a
    minusdifbplusa = - difb[0] + a

    #+a
    charge1 = copy.deepcopy(charge)
    for i in range(len(charge1)):
        charge1[i].coords[0] += a
        newcharge.append(charge_xyz(charge1[i].coords, charge1[i].charge))

    #-a
    charge2 = copy.deepcopy(charge)
    for i in range(len(charge2)):
        charge2[i].coords[0] -= a
        newcharge.append(charge_xyz(charge2[i].coords, charge2[i].charge))

    #+b
    charge3 = copy.deepcopy(charge)
    for i in range(len(charge3)):
         charge3[i].coords[0] += difb[0]
         charge3[i].coords[1] += difb[1]
         newcharge.append(charge_xyz(charge3[i].coords, charge3[i].charge))
    
    #+b+a
    charge4 = copy.deepcopy(charge)
    for i in range(len(charge4)): 
        charge4[i].coords[0] += difbplusa         
        charge4[i].coords[1] += difb[1] 
        newcharge.append(charge_xyz(charge4[i].coords, charge4[i].charge))
    
    #+b-a
    charge5 = copy.deepcopy(charge)
    for i in range(len(charge5)): 
        charge5[i].coords[0] += difbminusa    
        charge5[i].coords[1] += difb[1]
        newcharge.append(charge_xyz(charge5[i].coords, charge5[i].charge))

    #-b
    charge6 = copy.deepcopy(charge)
    for i in range(len(charge6)):
        charge6[i].coords[0] -= difb[0]
        charge6[i].coords[1] -= difb[1]
        newcharge.append(charge_xyz(charge6[i].coords, charge6[i].charge))

    #-b-a
    charge7 = copy.deepcopy(charge)
    for i in range(len(charge7)):
        charge7[i].coords[0] += minusdifbminusa
        charge7[i].coords[1] -= difb[1]
        newcharge.append(charge_xyz(charge7[i].coords, charge7[i].charge))

    #-b+a
    charge8 = copy.deepcopy(charge)
    for i in range(len(charge8)):
        charge8[i].coords[0] += minusdifbplusa
        charge8[i].coords[1] -= difb[1]
        newcharge.append(charge_xyz(charge8[i].coords, charge8[i].charge))

    #Begin second layer
    #+2a
    charge9 = copy.deepcopy(charge)
    for i in range(len(charge9)):
        charge9[i].coords[0] += 2*a
        newcharge.append(charge_xyz(charge9[i].coords, charge9[i].charge))

    #-2a
    charge10 = copy.deepcopy(charge)
    for i in range(len(charge10)):
        charge10[i].coords[0] -= 2*a
        newcharge.append(charge_xyz(charge10[i].coords, charge10[i].charge))

    #-2a+b
    charge11 = copy.deepcopy(charge)
    for i in range(len(charge11)):
        charge11[i].coords[0] += -2*a + difb[0]
        charge11[i].coords[1] += difb[1]
        newcharge.append(charge_xyz(charge11[i].coords, charge11[i].charge))

    #-2a+2b
    charge12 = copy.deepcopy(charge)
    for i in range(len(charge12)):
        charge12[i].coords[0] += -2 * (a + difb[0])
        charge12[i].coords[1] += 2*difb[1]
        newcharge.append(charge_xyz(charge12[i].coords, charge12[i].charge))

    #-a+2b
    charge13 = copy.deepcopy(charge)
    for i in range(len(charge13)):
        charge13[i].coords[0] += -a + 2*difb[0]
        charge13[i].coords[1] += 2*difb[1]
        newcharge.append(charge_xyz(charge13[i].coords, charge13[i].charge))

    #+a+2b
    charge14 = copy.deepcopy(charge)
    for i in range(len(charge14)):
        charge14[i].coords[0] += a + 2*difb[0]
        charge14[i].coords[1] += 2*difb[1]
        newcharge.append(charge_xyz(charge14[i].coords, charge14[i].charge))

    #+2a+2b
    charge15 = copy.deepcopy(charge)
    for i in range(len(charge15)):
        charge15[i].coords[0] += 2*(a + difb[0])
        charge15[i].coords[1] += 2*difb[1]
        newcharge.append(charge_xyz(charge15[i].coords, charge15[i].charge))

    #+2a+b
    charge16 = copy.deepcopy(charge)
    for i in range(len(charge16)):
        charge16[i].coords[0] += 2*a + difb[0]
        charge16[i].coords[1] += difb[1]
        newcharge.append(charge_xyz(charge16[i].coords, charge16[i].charge))

    #+2b
    charge17 = copy.deepcopy(charge)
    for i in range(len(charge17)):
        charge17[i].coords[0] += 2*difb[0]
        charge17[i].coords[1] += 2*difb[1]
        newcharge.append(charge_xyz(charge14[i].coords, charge14[i].charge))

    #+2b-a
    charge18 = copy.deepcopy(charge)
    for i in range(len(charge18)):
        charge14[i].coords[0] += -a + 2*difb[0]
        charge14[i].coords[1] += 2*difb[1]
        newcharge.append(charge_xyz(charge14[i].coords, charge14[i].charge))

    #-2a-b
    charge19 = copy.deepcopy(charge)
    for i in range(len(charge19)):
        charge19[i].coords[0] += -2*a - difb[0]
        charge19[i].coords[1] -= difb[1]
        newcharge.append(charge_xyz(charge19[i].coords, charge19[i].charge))

    #+2a-b
    charge20 = copy.deepcopy(charge)
    for i in range(len(charge20)):
        charge20[i].coords[0] += 2*a - difb[0]
        charge20[i].coords[1] -= difb[1]
        newcharge.append(charge_xyz(charge20[i].coords, charge20[i].charge))

    #+2a-2b
    charge21 = copy.deepcopy(charge)
    for i in range(len(charge21)):
        charge21[i].coords[0] += 2*(a - difb[0])
        charge21[i].coords[1] -= 2*difb[1]
        newcharge.append(charge_xyz(charge21[i].coords, charge21[i].charge))

    #+a-2b
    charge22 = copy.deepcopy(charge)
    for i in range(len(charge22)):
        charge22[i].coords[0] += a - 2*difb[0]
        charge22[i].coords[1] -= 2*difb[1]
        newcharge.append(charge_xyz(charge22[i].coords, charge22[i].charge))

    #-2b
    charge23 = copy.deepcopy(charge)
    for i in range(len(charge23)):
        charge23[i].coords[0] -= 2*difb[0]
        charge23[i].coords[1] -= 2*difb[1]
        newcharge.append(charge_xyz(charge23[i].coords, charge23[i].charge))

    #-2b-a
    charge24 = copy.deepcopy(charge)
    for i in range(len(charge24)):
        charge24[i].coords[0] -= 2*difb[0] + a
        charge24[i].coords[1] -= 2*difb[1]
        newcharge.append(charge_xyz(charge24[i].coords, charge24[i].charge))

    #-2b-2a
    charge25 = copy.deepcopy(charge)
    for i in range(len(charge25)):
        charge25[i].coords[0] -= 2*difb[0] + a
        charge25[i].coords[1] -= 2*difb[1]
        newcharge.append(charge_xyz(charge25[i].coords, charge25[i].charge))

    return newcharge

'''Funktion shiftcoord:
    verschiebt die eingegebene Koordinate um den eingegebenen Wert
    und speichert die neue Koordinate zusammen mit den unveränderten 
    in einer neuen Liste.'''
def shiftcoord(xyz, input_file, shiftsize):
    newcoord = []
    coord = []
    
    for i in range(len(input_file)):
        if xyz == "x":
            input_file[i].coords[0] = (input_file[i].coords[0] + shiftsize)
            newcoord.append(atom_xyz(input_file[i].symbol, input_file[i].coords))
              
        elif xyz == "y":
            input_file[i].coords[1] = (input_file[i].coords[1] + shiftsize)
            newcoord.append(atom_xyz(input_file[i].symbol, input_file[i].coords))
            
        elif xyz == "z":
            input_file[i].coords[2] = (input_file[i].coords[2] + shiftsize)
            newcoord.append(atom_xyz(input_file[i].symbol, input_file[i].coords)) 

    return newcoord

'''Funktion shiftcoordcharges:
    verschiebt die eingegebene Koordinate um den eingegebenen Wert
    und speichert die neue Koordinate zusammen mit den unveränderten 
    in einer neuen Liste.'''
def shiftcoordcharges(xyz,input_file,shiftsize):
    newcoord = []
    coord = []
    
    for i in range(len(input_file)):
        if xyz == "x":
            input_file[i].coords[0] = (input_file[i].coords[0] + shiftsize)
            newcoord.append(charge_xyz(input_file[i].coords,input_file[i].charge))
            
        elif xyz == "y":
            input_file[i].coords[1] = (input_file[i].coords[1] + shiftsize)
            newcoord.append(charge_xyz(input_file[i].coords, input_file[i].charge))
            
        elif xyz == "z":
            input_file[i].coords[2] = (input_file[i].coords[2] + shiftsize)
            newcoord.append(charge_xyz(input_file[i].coords, input_file[i].charge)) 
    
    return newcoord

            
'''Mathematische Umformung um in Richtung der b-Achse des PDIRs zu verschieben.
	Gibt eine Liste zurück mit x- und y-Koordinaten,
	welche die Differenz zur b-Achse beschreibt'''
def getbPDIR(gamma, b):
    difb = []

    gammastrich = 180 - gamma
    g = m.radians(gammastrich)

    y = m.sin(g) * b
    x = m.cos(g) * b

    difb.append(x)
    difb.append(y)

    return difb

'''Mathematische Umformung um in Richtung der c-Achse des PDIRs zu verschieben.
    Gibt eine Liste zurück mit x-, y- und z-Koordinaten,
    welche die Differenz zur c-Achse beschreibt'''
def getcPDIR(alpha, beta, gamma, c,):
    difc = []

    betastrich = 90 - beta
    gammastrich = 180 - gamma

    g = m.radians(gammastrich)
    betastra = m.radians(betastrich)
    alphara = m.radians(alpha)
    betara = m.radians(beta)

    x = m.cos(betara) * c
    y = m.cos(alphara) * c * m.sin(g)
    z = m.cos(betastra) * c

    difc.append(x)
    difc.append(y)
    difc.append(z)

    return difc

'''Funktion rot_axis_angle:
    Erstellt eine Rotationsmatrix'''
def rot_axis_angle(axis, angle): # [0,1,0] |90
    axis /= np.linalg.norm(axis)
    theta = m.radians(angle)
    costheta = m.cos(theta)
    sintheta = m.sin(theta)
    Ex = np.array(
        [[0.0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0.0]])
    A = costheta * \
        np.eye(3) + (1.0 - costheta) * \
        np.dot(axis.reshape(3, 1), axis.reshape(1, 3)) + \
        sintheta * Ex
    return A

'''Funktion product_matrix_vector:
    Multipliziert eine Matrix mit einem Vektor und gibt
    das Ergebnis zurück.
    Speichert das Ergebnis in der Klasse atom.'''
def product_matrix_vector(matrix, input_file):
    newcoord = []

    for i in range(len(input_file)):
        vec = input_file[i].coords
        newvec = np.dot(matrix,vec)
        newcoord.append(atom_xyz(input_file[i].symbol, newvec))

    return newcoord

'''Funktion product_matrix_vector_charges:
    Multipliziert eine Matrix mit einem Vektor und gibt das Ergebnis zurück.
    Speichert das Ergebnis in der Klasse charge.'''
def product_matrix_vector_charges(matrix, input_file):
    newcoord = []

    for i in range(len(input_file)):
        vec = input_file[i].coords
        newvec = np.dot(matrix,vec)
        newcoord.append(charge_xyz(newvec, input_file[i].charge))

    return newcoord

'''Speichert die Symbole und Koordinaten von DIP und PDIR in einer Liste ab.'''
def getonefile (retPDIR, ret):
    ret1 = copy.deepcopy(ret)

    for i in range(len(retPDIR)):
        ret1.append(atom_xyz(retPDIR[i].symbol, retPDIR[i].coords))

    return ret1

'''Schreibt sowohl die Ladungen (mit Atomsymbol) als auch
	DIP und PDIR in eine xyz Datei.'''
def outputeverythingDIPPDIR(retDIP, retPDIR, retDIPPDIR, chargesDIP,
    chargesPDIR, step, step2):
    totallength = int(len(retDIPPDIR))
    lengthchargesDIP = int(len(chargesDIP))
    lengthchargesPDIR = int(len(chargesPDIR))
    ret1 = copy.deepcopy(retDIP)
    ret2 = []
    ret3 = []

    lengthfile = lengthchargesDIP + totallength + lengthchargesPDIR

    f = open("everything{}_{}.xyz".format(step, step2), "w")
    f.write(str(lengthfile))
    f.write("\n")
    f.write("\n")

    i = 0
    while i < 6:
        for j in range(len(ret1)):
            ret2.append(atom_xyz(ret1[j].symbol, ret1[j].coords))
        i += 1

    n = 0
    while n < 8:
        for j in range(len(retPDIR)):
            ret3.append(atom_xyz(retPDIR[j].symbol, retPDIR[j].coords))
        n += 1

    for j in range(len(chargesDIP)):
        f.write(str(ret2[j].symbol))
        f.write("    ")
        k = 0 
        while k < 3:
            f.write(str(np.array(chargesDIP[j].coords[k])))
            f.write("   ")
            k += 1
        f.write("\n")

    for j in range(len(chargesPDIR)):
        f.write(str(ret3[j].symbol))
        f.write("    ")
        k = 0 
        while k < 3:
            f.write(str(np.array(chargesPDIR[j].coords[k])))
            f.write("   ")
            k += 1
        f.write("\n")
    
    for l in range(len(retDIPPDIR)):
        f.write(str(retDIPPDIR[l].symbol))
        f.write("   ")
        m = 0
        while m < 3:
            f.write(str(np.array(retDIPPDIR[l].coords[m])))
            f.write("   ")
            m += 1
        f.write("\n")

    f.close()


'''Funktion output:
    Öffnet eine neue xyz-Datei und speichert dort
    die veränderten Koordinaten und Symbole ab.'''
def output(newcoord, step, step2):
    f = open("neueKoordinaten{}_{}.xyz".format(step, step2), "w")
    f.write(str(int(len(newcoord))))
    f.write("\n")
    f.write("\n")
    for a in range(len(newcoord)):
        f.write(str(newcoord[a].symbol))
        f.write("    ")
        j = 0
        while j < 3:
            f.write(str(np.array(newcoord[a].coords[j])))
            f.write("   ")
            j += 1
        f.write("\n")
    f.close()

'''Öffnet eine neue xyz-Datei und speichert dort
    die veränderten Koordinaten und Ladungen ab.'''
def outputchargesDIP(newcoord, step, step2):
    f = open("neuechargesDIP{}_{}.txt".format(step, step2), "w")
    f.write(str(int(len(newcoord))))
    f.write("\n")
    f.write("\n")
    for a in range(len(newcoord)):
        j = 0
        while j<3:
            f.write(str(np.array(newcoord[a].coords[j])))
            f.write("   ")
            j += 1
        f.write(str(float(newcoord[a].charge)))
        f.write("\n")
    f.close()

'''Funktion make_com:
    Erstellt eine com-Datei'''
def make_com(xyz1, xyz2, chargesDIP, chargesPDIR, char):
    if char == "Ja":
        head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
#p wB97XD/cc-pVDZ charge gfinput pop=full TD(nstates=8) 

Dimer

0 1
"""
        with open("gauss.com","w") as out_file:
            out_file.write(head)
            for struct in xyz1:
                out_file.write("{} {:12.6f} {:12.6f} {:12.6f}\n".format(struct.symbol, *struct.coords))
            for struct in xyz2:
                out_file.write("{} {:12.6f} {:12.6f} {:12.6f}\n".format(struct.symbol, *struct.coords))
            out_file.write("\n")
            for struct in chargesDIP:
                out_file.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(*struct.coords, struct.charge))
            for struct in chargesPDIR:
                out_file.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(*struct.coords, struct.charge))
            out_file.write("\n\n\n\n")

    else:
        head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
#p wB97XD/cc-pVDZ gfinput pop=full TD(nstates=8) 

Dimer

0 1
"""
        with open("gauss.com","w") as out_file:
            out_file.write(head)
            for struct in xyz1:
                out_file.write("{} {:12.6f} {:12.6f} {:12.6f}\n".format(struct.symbol, *struct.coords))
            for struct in xyz2:
                out_file.write("{} {:12.6f} {:12.6f} {:12.6f}\n".format(struct.symbol, *struct.coords))
            out_file.write("\n")
    

'''Funktion make_sh:
    Erstellt eine sh-Datei'''
def make_sh(step, step2, xyz1, xyz2):
    head = """#PBS -l nodes=1:ppn=8
#PBS -l walltime=800:00:00
#PBS -l mem=8GB
#PBS -N {}_{}_{}_{}Dimer

cd $PBS_O_WORKDIR
cp -rf * $TMPDIR
cd $TMPDIR

# setenv g16root /apps/gaussian/g16_A
# source $g16root/g16/bsd/g16.login

# setenv GAUSS_SCRDIR $TMPDIR

g16  gauss.com
rm Gau-*

cp -rf * $PBS_O_WORKDIR

"""
    with open("gauss.sh","w") as out_file:
        out_file.write(head.format("charge-field", step, step2, xyz1, xyz2))

#Funktion die alle Usereingaben in eine Datei schreibt.
def useroutput(char, DIP2, dup, dup2, verschiebung1, verschiebung2,
    xyz1, xyz2, shiftstart, shiftstart2, shiftsize, shiftsize2, shiftlength,
    shiftlength2, wahl, axis, angle):
    f = open("Usereingaben.txt","w")
    f.write("Folgende Usereingaben wurden bei dieser Rechnung benutzt: ")
    f.write("\n\n")
    f.write("Wurden Ladungen berücksichtigt: %s" %char)
    f.write("\n")
    f.write("\n")
    f.write("Ist das 2. DIP-Molekül aus der Einheitszelle dabei: %s" %DIP2)
    f.write("\n")
    f.write("\n")
    f.write("Wurde DIP dupliziert: %s" %dup)
    f.write("\n")
    f.write("Um welche Koordinate wurde es verschoben: %s" %verschiebung1)
    f.write("\n")
    f.write("\n")
    f.write("Wurde PDIR dupliziert: %s" %dup2)
    f.write("\n")
    f.write("Um welche Koordinate wurde es verschoben: %s" %verschiebung2)
    f.write("\n")
    f.write("\n")
    f.write("Die DIP-Moleküle wurden um die Koordinate %s verschoben." %xyz1)
    f.write("\n")
    f.write("Die DIP-Moleküle wurden um die Koordinate %s verschoben." %xyz2)
    f.write("\n")
    f.write("1. Koordinate:")
    f.write("\n")
    f.write("   Startpunkt: %s" %shiftstart)
    f.write("\n")
    f.write("   Schrittlänge: %s" %shiftsize)
    f.write("\n")
    f.write("   Gesamtlänge: %s" %shiftlength)
    f.write("\n")
    f.write("2. Koordinate:")
    f.write("\n")
    f.write("   Startpunkt: %s" %shiftstart2)
    f.write("\n")
    f.write("   Schrittlänge: %s" %shiftsize2)
    f.write("\n")
    f.write("   Gesamtlänge: %s" %shiftlength2)
    f.write("\n")
    f.write("\n")
    f.write("Wurde der DIP-Kristall rotiert: %s" %wahl)
    f.write("\n")
    f.write("   Achse: %s" %axis)
    f.write("\n")
    f.write("   Winkel: %s" %angle)
    f.write("\n")
    f.close()

def thirdAxis(first, second):
    if first == "x":
        if second == "y":
            return "z"
        if second == "z":
            return "y"
    elif first == "y":
        if second == "x":
            return "z"
        if second == "z":
            return "x"
    elif first == "z":
        if second == "x":
            return "y"
        if second == "y":
            return x






'''Funktion main:
    Beinhaltet alle Ausführungen des Codes.'''
def main():
    #Initalisieren der Variablen
    chargesDIP = []
    new_charges = []
    first_charges = []
    geo_chargesPDIR = []

    first_input = []
    verschiebung1 = ""
    verschiebung2 = ""
    axis = ""
    angle = ""
    step = 0
    step2 = 0

    #Länge der a-, b- und c-Achse des DIP-Einheitskristalles
    lengthaDIP = 7.1709
    lengthbDIP = 8.5496
    lengthcDIP = 16.7981

    #Länge der a-, b- und c-Achse des PDIR-Einheitskristalles
    lengthaPDIR = 6.3022
    lengthbPDIR = 8.80954
    lengthcPDIR = 16.5571
    
    #Winkel des PDIR-Einheitskristalls
    alphaPDIR = 88.3219
    betaPDIR = 84.604
    gammaPDIR = 110.631
    
    #Differenz zwischen y- und b-Achse und z- und c-Achse
    difb = getbPDIR(gammaPDIR, lengthbPDIR)
    difc = getcPDIR(alphaPDIR, betaPDIR, gammaPDIR, lengthbPDIR)

    #Einlesen der Dateien
    file_input = readXYZ("output_first_DIP_molecule___with_opt_geo.xyz")
    file_input1 = readXYZ("output_second_DIP_molecule___with_opt_geo.xyz")
    file_input2 = readXYZ("PDIR_S0_1.xyz")

    #Initalisieren der Variablen mit Userinput
    print("VORSICHT: Das Molekül wird in das Center of Geometry geschoben!")

    char = str(input("Sollen Ladungen berücksichtigt werden? (Ja, Nein): "))
    DIP2 = input("Soll das 2. DIP-Molekül aus der Einheitszelle dazugeladen werden? (Ja, Nein): ")
    dup = input("Soll DIP dupliziert werden? (Ja, Nein): ")
    dup2 = input("Soll PDIR dupliziert werden? (Ja, Nein): ")
    
    #Liest das 2. DIP-Molekül ein und hängt es an das erste an.
    if DIP2 == "Ja":
        for i in range(len(file_input)):
            file_input.append(atom_xyz(file_input1[i].symbol, file_input1[i].coords))

    '''Liest die Ladungen von DIP und PDIR ein, verschiebt sie,
    	speichert sie in 2 Listen und macht centerofGeo von beiden'''
    file_geo = getcenterofGeo(file_input)
    if char == "Ja":
        chargesDIP1 = readcharges("first_DIP_charges.txt")
        chargesDIP2 = readcharges("second_DIP_charges.txt")
        chargesPDIR = readcharges("PDIR_S0_1_charges.txt")

        geo_chargesPDIR = getcenterofGeo(chargesPDIR)
    
        shiftchargesdip1 = getchargesDIP1(chargesDIP1, lengthaDIP, lengthbDIP) 
        shiftchargesdip2 = getchargesDIP2(chargesDIP2, lengthaDIP, lengthbDIP)
        PDIRcharges = getchargesPDIR(chargesPDIR, difb, lengthaPDIR)

        for i in range(len(shiftchargesdip2)):
            shiftchargesdip1.append(charge_xyz(shiftchargesdip2[i].coords, shiftchargesdip2[i].charge))
        
        geo_chargesDIP = getcenterofGeo(shiftchargesdip1)
        geo_chargesPDIR = PDIRcharges



        
        
      
    #Erstellt Center of Geometry der Moleküle
    
    file_geo1 = getcenterofGeo(file_input1)
    file_geo2 = getcenterofGeo(file_input2)  

    #Duplizieren sowohl von DIP als auch von PDIR
    if dup == "Ja":
        verschiebung1 = input("In welche Richtung soll das zweite DIP-Molekül verschoben werden? (xp, xn, yp, yn, zp, zn):")
        file_geo_verschoben = duplicateDIP(file_geo, verschiebung1, lengthaDIP, lengthbDIP)
        for i in range(len(file_geo_verschoben)):
            file_geo.append(atom_xyz(file_geo_verschoben[i].symbol, file_geo_verschoben[i].coords))

    if dup2 == "Ja":
        verschiebung2 = input("In welche Richtung soll das zweite PDIR-Molekül verschoben werden? (ap, an, bp, bn, cp, cn):") 
        file_geo2_verschoben = duplicatePDIR(file_geo2, lengthaPDIR, difb, difc, verschiebung2)
        for i in range(len(file_geo2_verschoben)):
            file_geo2.append(atom_xyz(file_geo2_verschoben[i].symbol, file_geo2_verschoben[i].coords))

    #Rotation des Kristalls
    wahl = input("Willst du den Kristall rotieren? (Ja, Nein): ")
    if wahl == "Ja":
        axis = input("Um welche Achse soll rotiert werden? (a, b, c): ")
        if axis == "a":
            axis = [1,0,0]
            angle = 180 - alphaPDIR
            print("Vorsicht: y = z und z = y")
        elif axis == "b":
            axis = [0,1,0]
            angle = 180 - betaPDIR
            print("Vorsicht: x = z und z = x")
        elif axis == "c":
            axis = [0,0,1]
            angle = 180 - gammaPDIR
            print("Vorsicht: x = y und y = x")

        file_geo = product_matrix_vector(rot_axis_angle(axis,angle),file_geo)

        if char == "Ja":
            geo_chargesDIP = product_matrix_vector_charges(rot_axis_angle(axis,angle),geo_chargesDIP)

    #Initalisieren weiterer Varibalen
    xyz1 = input("Gib die 1. Koordinate an, die verschoben werden soll (x, y, z): ")
    xyz2 = input("Gib die 2. Koordinate an, die verschoben werden soll (x, y, z): ")
    distance = float(input("Gib den Abstand für die dritte Koordinate an: "))

    shiftstart = float(input("Gib den Startpunkt der ersten Koordinate an: "))
    shiftsize = float(input("Gib die Schrittweite der Verschiebung der ersten Koordinate an: "))
    shiftlength = float(input("Gib den Endpunkt der erste Koordinate an: "))
    shiftstart2 = float(input("Gib den Startpunkt der zweiten Koordinate an: "))
    shiftsize2 = float(input("Gib die Schrittweite der Verschiebung der zweiten Koordinate an: "))
    shiftlength2 = float(input("Gib den Endpunkt der zweiten Koordinate an: "))


    #Erstellt einen neuen Ordner, geht hinein und macht ihn zum Hauptordner
    new= input("Ordnername: ")
    os.mkdir(new)
    os.chdir(new)
    main_dir = os.getcwd()

    #Erstellt eine Datei, welche die Usereingaben beinhalten.
    useroutput(char, DIP2, dup, dup2, verschiebung1, verschiebung2,
        xyz1, xyz2, shiftstart, shiftstart2, shiftsize, shiftsize2,
        shiftlength, shiftlength2, wahl, axis, angle)

    #Erstellt einen neuen Ordner und geht hinein
    newdir="{}_{}_vacuum".format(step, xyz1)
    os.mkdir(newdir)
    os.chdir(newdir)

    #Erstellt einen neuen Ordner und geht hinein
    newdir2 = "{}_{}_vacuum".format(step2, xyz2)
    os.mkdir(newdir2)
    os.chdir(newdir2)

    '''Berechnet den Startpunkt, speichert die Symbole und
    	veränderten Koordinaten in eine Datei und liest sie wieder ein.'''
    shiftxyz1 = shiftcoord(xyz1, file_geo, shiftstart)
    shiftz = shiftcoord(thirdAxis(xyz1, xyz2), shiftxyz1, distance)
    shiftall = shiftcoord(xyz2, shiftz, shiftstart2)

    output(shiftall, step, step2)

    first_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))

    '''Wenn Ladungen berücksichtigt werden:
    	Verändert die Koordinaten der Ladungen zum Startpunkt, speichert sie in einer Datei
    	liest sie ein und gibt eine xyz-Datei aus, in der die Ladungen und die beiden
    	Moleküle DIP und PDIR vorhanden sind.'''
    if char == "Ja":
        shiftchargexyz1 = shiftcoordcharges(xyz1, geo_chargesDIP, shiftstart)
        shiftchargez = shiftcoordcharges(thirdAxis(xyz1, xyz2), shiftchargexyz1, distance)
        shiftallcharges = shiftcoordcharges(xyz2, shiftchargez, shiftstart2)

        outputchargesDIP(shiftallcharges, step, step2)

        first_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
        
        DIPPDIR = getonefile(file_geo2, first_input)

        outputeverythingDIPPDIR(first_input, file_geo2, DIPPDIR, first_charges, geo_chargesPDIR, step, step2)

    #Erstellen der com- und sh-Datei
    make_com(first_input, file_geo2, first_charges, geo_chargesPDIR, char)
    make_sh(step, step2, xyz1, xyz2)

    #Abschicken der Rechnung
    subprocess.call('qsub -V gauss.sh', shell = True)

    #Schleife, die die ebengenannten Schritte wiederholt.
    i = shiftstart
    while i <= shiftlength:

        j = shiftstart2
        while j < shiftlength2:

            new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))

            if char == "Ja":
                new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
                
            step2 = step2 + 1

            os.chdir(main_dir)
            os.chdir(newdir)

            newdir3 = "{}_{}_vacuum".format(step2, xyz2)
            os.mkdir(newdir3)
            os.chdir(newdir3)

            shifty = shiftcoord(xyz2, new_input, shiftsize2)

            output(shifty, step, step2)

            if char == "Ja":
                shiftchargey = shiftcoordcharges(xyz2, new_charges, shiftsize2)

                outputchargesDIP(shiftchargey, step, step2)

                new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))
                new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
                DIPPDIR = getonefile(file_geo2, new_input)

                outputeverythingDIPPDIR(new_input, file_geo2, DIPPDIR, new_charges, geo_chargesPDIR, step, step2)
            
 
            make_com(new_input, file_geo2, new_charges, geo_chargesPDIR, char)
            make_sh(step, step2, xyz1, xyz2)

            subprocess.call('qsub -V gauss.sh', shell = True)

            j = j + shiftsize2
        
        if i < shiftlength:
            step2 = 0

            os.chdir(main_dir)
            os.chdir(newdir)
            os.chdir(newdir2)

            new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))

            if char == "Ja":
                new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))
                new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))

            step = step + 1
        
            os.chdir(main_dir)

            newdir = "{}_{}_vacuum".format(step, xyz1)
            os.mkdir(newdir)
            os.chdir(newdir)

            newdir5 = "{}_{}_vacuum".format(step2, xyz2)
            os.mkdir(newdir5)
            os.chdir(newdir5)

            shiftx = shiftcoord(xyz1, new_input, shiftsize)

            output(shiftx, step, step2)

            if char == "Ja":
                shiftchargex = shiftcoordcharges(xyz1, new_charges, shiftsize)

                outputchargesDIP(shiftchargex, step, step2)

                new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))
                new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
                DIPPDIR = getonefile(file_geo2, new_input)
                
                outputeverythingDIPPDIR(new_input, file_geo2, DIPPDIR, new_charges, geo_chargesPDIR, step, step2)

            make_com(new_input, file_geo2, new_charges, geo_chargesPDIR, char)
            make_sh(step, step2, xyz1, xyz2)

            subprocess.call('qsub -V gauss.sh', shell = True)

        else:
        	break
   
        i = i+shiftsize
    
if __name__=="__main__": main()   





    
