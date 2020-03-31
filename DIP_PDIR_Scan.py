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

'''Funktion moveToCenterofGeo:
    Berechnet Center of Geometry und gibt die
    veränderten Koordinaten zurück.'''

def moveToCenterofGeo(ret):

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

def getCenterOfGeo(ret):
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
    
    return centerOfGeo


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
countChargesDIP1 = 0
def getchargesDIP1(charge, a, b, c, numberLayers, axisDIP, xyz, reduceChargeSize):
    global countChargesDIP1
    newcharge =[]
    
    shiftb = True # Should charges be duplicated along b axis (before rotating)
    shifta = True
    if axisDIP == "a" and (xyz == "z" or xyz == "x"):
        shiftb = False
    if axisDIP == "b" and (xyz == "z" or xyz == "y"):
        shifta = False

    full = True
    if reduceChargeSize == "Ja":
        full = False
    
    # First layer DIP1
    if numberLayers >= 1:
        if shifta:
            countChargesDIP1 += 1
            charge1 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge1[i].coords[0] += a
                newcharge.append(charge_xyz(charge1[i].coords, charge1[i].charge))

        if shifta:
            countChargesDIP1 += 1
            charge2 = copy.deepcopy(charge)
            for i in range(len(charge)):   
                charge2[i].coords[0] -= a
                newcharge.append(charge_xyz(charge2[i].coords, charge2[i].charge))

        if shifta and full:
            countChargesDIP1 += 1
            charge3 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge3[i].coords[0] += 2*a
                newcharge.append(charge_xyz(charge3[i].coords, charge3[i].charge))

        if shifta and full:
            countChargesDIP1 += 1
            charge4 = copy.deepcopy(charge)
            for i in range(len(charge)):   
                charge4[i].coords[0] -= 2*a
                newcharge.append(charge_xyz(charge4[i].coords, charge4[i].charge))

        # Begin row DIP1 +b
        if shiftb:
            countChargesDIP1 += 1
            charge5 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge5[i].coords[1] += b
                newcharge.append(charge_xyz(charge5[i].coords, charge5[i].charge))
        
        if shiftb and shifta:
            countChargesDIP1 += 1
            charge6 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge6[i].coords[0] += a
                charge6[i].coords[1] += b
                newcharge.append(charge_xyz(charge6[i].coords, charge6[i].charge))

        if shiftb and shifta:
            countChargesDIP1 += 1
            charge7 = copy.deepcopy(charge)
            for i in range(len(charge)):   
                charge7[i].coords[0] -= a
                charge7[i].coords[1] += b
                newcharge.append(charge_xyz(charge7[i].coords, charge7[i].charge))

        if shiftb and shifta and full:
            countChargesDIP1 += 1
            charge8 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge8[i].coords[0] += 2*a
                charge8[i].coords[1] += b
                newcharge.append(charge_xyz(charge8[i].coords, charge8[i].charge))

        if shiftb and shifta and full:
            countChargesDIP1 += 1
            charge9 = copy.deepcopy(charge)
            for i in range(len(charge)):   
                charge9[i].coords[0] -= 2*a
                charge9[i].coords[1] += b
                newcharge.append(charge_xyz(charge9[i].coords, charge9[i].charge))

        # Begin row DIP1 -b
        if shiftb:
            countChargesDIP1 += 1
            charge10 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge10[i].coords[1] -= b
                newcharge.append(charge_xyz(charge10[i].coords, charge10[i].charge))
        
        if shiftb and shifta:
            countChargesDIP1 += 1
            charge11 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge11[i].coords[0] += a
                charge11[i].coords[1] -= b
                newcharge.append(charge_xyz(charge11[i].coords, charge11[i].charge))

        if shiftb and shifta:
            countChargesDIP1 += 1
            charge12 = copy.deepcopy(charge)
            for i in range(len(charge)):   
                charge12[i].coords[0] -= a
                charge12[i].coords[1] -= b
                newcharge.append(charge_xyz(charge12[i].coords, charge12[i].charge))

        if shiftb and shifta and full:
            countChargesDIP1 += 1
            charge13 = copy.deepcopy(charge)
            for i in range(len(charge)):
                charge13[i].coords[0] += 2*a
                charge13[i].coords[1] -= b
                newcharge.append(charge_xyz(charge13[i].coords, charge13[i].charge))

        if shiftb and shifta and full:
            countChargesDIP1 += 1
            charge14 = copy.deepcopy(charge)
            for i in range(len(charge)):   
                charge14[i].coords[0] -= 2*a
                charge14[i].coords[1] -= b
                newcharge.append(charge_xyz(charge14[i].coords, charge14[i].charge))

    return newcharge

'''Bildet die Ladungen um DIP2 und speichert sie in einer Liste.'''
countChargesDIP2 = 0 
def getchargesDIP2(charge, a, b, c, numberLayers, DIP2, axisDIP, xyz, reduceChargeSize):
    global countChargesDIP2

    chargegeaendert = []
    duplicateA = True
    if axisDIP == "b" and (xyz == "z" or xyz == "y"):
        duplicateA = False

    full = True
    if reduceChargeSize == "Ja":
        full = False

    if not(axisDIP == "a" and (xyz == "z" or xyz == "x")):

        if DIP2 == "Nein" and duplicateA:
            countChargesDIP2 += 1
            for i in range(len(charge)):
                chargegeaendert.append(charge_xyz(charge[i].coords, charge[i].charge))

        if duplicateA and full:
            countChargesDIP2 +=1  
            charge1 = copy.deepcopy(charge)
            for i in range(len(charge1)):
                charge1[i].coords[0] -= 2 * a
                chargegeaendert.append(charge_xyz(charge1[i].coords, charge1[i].charge))

        countChargesDIP2 +=1
        charge2 = copy.deepcopy(charge)
        for i in range(len(charge2)):
            charge2[i].coords[0] -= a
            chargegeaendert.append(charge_xyz(charge2[i].coords, charge2[i].charge))
        
        if duplicateA:
            countChargesDIP2 +=1
            charge3 = copy.deepcopy(charge)
            for i in range(len(charge3)):
                charge3[i].coords[1] += b
                chargegeaendert.append(charge_xyz(charge3[i].coords, charge3[i].charge))

        if duplicateA and full:
            countChargesDIP2 +=1
            charge4 = copy.deepcopy(charge)
            for i in range(len(charge4)):
                charge4[i].coords[0] += a
                charge4[i].coords[1] += b
                chargegeaendert.append(charge_xyz(charge4[i].coords, charge4[i].charge))
        
        if duplicateA and full:
            countChargesDIP2 +=1
            charge5 = copy.deepcopy(charge)
            for i in range(len(charge5)):
                charge5[i].coords[1] += b
                charge5[i].coords[0] -= 2 * a
                chargegeaendert.append(charge_xyz(charge5[i].coords, charge5[i].charge))

        countChargesDIP2 +=1
        charge6 = copy.deepcopy(charge)
        for i in range(len(charge6)):
            charge6[i].coords[1] += b
            charge6[i].coords[0] -= a
            chargegeaendert.append(charge_xyz(charge6[i].coords, charge6[i].charge))

        if duplicateA and full:
            countChargesDIP2 +=1
            charge7 = copy.deepcopy(charge)
            for i in range(len(charge7)):
                charge7[i].coords[0] += a
                chargegeaendert.append(charge_xyz(charge7[i].coords, charge7[i].charge))

        # Row DIP2 +2b
        if duplicateA and full:
            countChargesDIP2 +=1
            charge8 = copy.deepcopy(charge)
            for i in range(len(charge8)):
                charge8[i].coords[1] += 2*b
                chargegeaendert.append(charge_xyz(charge8[i].coords, charge8[i].charge))

        if duplicateA and full:
            countChargesDIP2 +=1
            charge9 = copy.deepcopy(charge)
            for i in range(len(charge9)):
                charge9[i].coords[0] += a
                charge9[i].coords[1] += 2*b
                chargegeaendert.append(charge_xyz(charge9[i].coords, charge9[i].charge))

        if duplicateA and full:
            countChargesDIP2 +=1
            charge10 = copy.deepcopy(charge)
            for i in range(len(charge10)):
                charge10[i].coords[1] += 2*b
                charge10[i].coords[0] -= 2 * a
                chargegeaendert.append(charge_xyz(charge10[i].coords, charge10[i].charge))
        if full:
            countChargesDIP2 +=1
            charge11 = copy.deepcopy(charge)
            for i in range(len(charge11)):
                charge11[i].coords[1] += 2*b
                charge11[i].coords[0] -= a
                chargegeaendert.append(charge_xyz(charge11[i].coords, charge11[i].charge))

        # Row DIP2 -b
        if duplicateA and full:
            countChargesDIP2 +=1
            charge12 = copy.deepcopy(charge)
            for i in range(len(charge12)):
                charge12[i].coords[1] -= b
                chargegeaendert.append(charge_xyz(charge12[i].coords, charge12[i].charge))

        if duplicateA and full:
            countChargesDIP2 +=1
            charge13 = copy.deepcopy(charge)
            for i in range(len(charge13)):
                charge13[i].coords[0] += a
                charge13[i].coords[1] -= b
                chargegeaendert.append(charge_xyz(charge13[i].coords, charge13[i].charge))

        if duplicateA and full:
            countChargesDIP2 +=1
            charge14 = copy.deepcopy(charge)
            for i in range(len(charge14)):
                charge14[i].coords[1] -= b
                charge14[i].coords[0] -= 2 * a
                chargegeaendert.append(charge_xyz(charge14[i].coords, charge14[i].charge))

        if full:
            countChargesDIP2 +=1
            charge15 = copy.deepcopy(charge)
            for i in range(len(charge15)):
                charge15[i].coords[1] -= b
                charge15[i].coords[0] -= a
                chargegeaendert.append(charge_xyz(charge15[i].coords, charge15[i].charge))


    return chargegeaendert

'''Bildet die Ladungen um PDIR und speichert sie in einer Liste.
Die Ladungen werden anhand der Rotationsaxe ausgewählt.'''
countChargesPDIR=0
def getchargesPDIR(charge, difb, a, difc, rotAxisDIP, numberLayers, direction, xyz, reduceChargeSize):
    newcharge = []
    global countChargesPDIR
    difbplusa = difb[0] + a
    difbminusa = difb[0] - a
    minusdifbminusa = - difb[0] - a
    minusdifbplusa = - difb[0] + a

    beginMultiplierA = -2
    beginMultiplierB = -2
    moleculesA = 5
    moleculesB = 5

    # Set values that rotation and translation are correct
    if rotAxisDIP == "a" and (xyz == "z" or xyz == "x"):
        rotAxisDIP = "0"

    if rotAxisDIP == "b" and (xyz == "z" or xyz == "y"):
        rotAxisDIP = "0"


    # Adjust values
    if rotAxisDIP == "a":
        beginMultiplierB = 1 - numberLayers
        moleculesB = numberLayers

    if rotAxisDIP == "b":
        beginMultiplierA = 1 - numberLayers
        moleculesA = numberLayers

    if reduceChargeSize == "Ja" and rotAxisDIP == "0":
        beginMultiplierA += 1
        beginMultiplierB += 1
        moleculesA -= 2
        moleculesB -= 2

    if reduceChargeSize == "Ja" and rotAxisDIP == "a":
        beginMultiplierA += 1
        moleculesA -= 2

    if reduceChargeSize == "Ja" and rotAxisDIP == "b":
        beginMultiplierB += 1
        moleculesB -= 2

    #Begin first layer (middle)
    if (numberLayers >= 1 and (rotAxisDIP == "c" or rotAxisDIP == "0")) or (rotAxisDIP == "a") or (rotAxisDIP == "b"):
        for i in range(moleculesB):
            for j in range(moleculesA):
                # Check which charge should be left out:
                if not ((direction == "bn" and beginMultiplierA == -j and (beginMultiplierB + i) == -1) or
                        (direction == "bp" and beginMultiplierA == -j and (beginMultiplierB + i) == 1) or
                        (direction == "an" and (beginMultiplierA + j) == -1 and beginMultiplierB == -i) or
                        (direction == "ap" and (beginMultiplierA + j) == 1 and beginMultiplierB == -i)):
                    if (beginMultiplierB != -i or beginMultiplierA != -j):
                        countChargesPDIR +=1
                        chargeFirstLayer = copy.deepcopy(charge)
                        for k in range(len(chargeFirstLayer)):
                            chargeFirstLayer[k].coords[0] += beginMultiplierA*a + beginMultiplierB*difb[0] + j*a + i*difb[0]
                            chargeFirstLayer[k].coords[1] += beginMultiplierB*difb[1] + i*difb[1]
                            newcharge.append(charge_xyz(chargeFirstLayer[k].coords, chargeFirstLayer[k].charge))

    #Begin second layer (-z)
    if (numberLayers >= 2 and (rotAxisDIP == "c" or rotAxisDIP == "0")) or (rotAxisDIP == "a") or (rotAxisDIP == "b"):
        for i in range(moleculesB):
            for j in range(moleculesA):
                if not (direction == "cn" and beginMultiplierA == -j and beginMultiplierB == -i):
                    countChargesPDIR +=1
                    chargeSecondLayer = copy.deepcopy(charge)
                    for k in range(len(chargeSecondLayer)):
                        chargeSecondLayer[k].coords[0] += beginMultiplierA*a + beginMultiplierB*difb[0] + j*a + i*difb[0] - difc[0]
                        chargeSecondLayer[k].coords[1] += beginMultiplierB*difb[1] + i*difb[1] - difc[1]
                        chargeSecondLayer[k].coords[2] += -difc[2]
                        newcharge.append(charge_xyz(chargeSecondLayer[k].coords, chargeSecondLayer[k].charge))

    #Begin third layer (+z)
    if (numberLayers >= 3 and (rotAxisDIP == "c" or rotAxisDIP == "0")) or (rotAxisDIP == "a") or (rotAxisDIP == "b"):
        for i in range(moleculesB):
            for j in range(moleculesA):
                if not (direction == "cp" and beginMultiplierA == -j and beginMultiplierB == -i):
                    countChargesPDIR +=1
                    chargeThirdLayer = copy.deepcopy(charge)
                    for k in range(len(chargeThirdLayer)):
                        chargeThirdLayer[k].coords[0] += beginMultiplierA*a + beginMultiplierB*difb[0] + j*a + i*difb[0] + difc[0]
                        chargeThirdLayer[k].coords[1] += beginMultiplierB*difb[1] + i*difb[1] + difc[1]
                        chargeThirdLayer[k].coords[2] += difc[2]
                        newcharge.append(charge_xyz(chargeThirdLayer[k].coords, chargeThirdLayer[k].charge))

    

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
def shiftcoordcharges(xyz, input_file, shiftsize):
    newcoord = []
    coord = []
    
    for i in range(len(input_file)):
        if xyz == "x":
            input_file[i].coords[0] = (input_file[i].coords[0] + shiftsize)
            newcoord.append(charge_xyz(input_file[i].coords, input_file[i].charge))
            
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
def getcPDIR(alpha, beta, gamma, c):
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
def outputeverythingDIPPDIR1D(retDIP, retPDIR, retDIPPDIR, chargesDIP,
    chargesPDIR, step):
    totallength = int(len(retDIPPDIR))
    lengthchargesDIP = int(len(chargesDIP))
    lengthchargesPDIR = int(len(chargesPDIR))
    ret1 = copy.deepcopy(retDIP)
    ret2 = []
    ret3 = []
    global countChargesPDIR
    global countChargesDIP1
    global countChargesDIP2

    lengthfile = lengthchargesDIP + totallength + lengthchargesPDIR

    f = open("everything{}.xyz".format(step), "w")
    f.write(str(lengthfile))
    f.write("\n")
    f.write("\n")

    i = 0
    while i < (countChargesDIP2+countChargesDIP1):
        for j in range(len(ret1)):
            ret2.append(atom_xyz(ret1[j].symbol, ret1[j].coords))
        i += 1

    n = 0
    while n < countChargesPDIR:
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

def outputeverythingDIPPDIR2D(retDIP, retPDIR, retDIPPDIR, chargesDIP,
    chargesPDIR, step, step2):
    totallength = int(len(retDIPPDIR))
    lengthchargesDIP = int(len(chargesDIP))
    lengthchargesPDIR = int(len(chargesPDIR))
    ret1 = copy.deepcopy(retDIP)
    ret2 = []
    ret3 = []
    global countChargesPDIR
    global countChargesDIP1
    global countChargesDIP2

    lengthfile = lengthchargesDIP + totallength + lengthchargesPDIR

    f = open("everything{}.xyz".format(step), "w")
    f.write(str(lengthfile))
    f.write("\n")
    f.write("\n")

    i = 0
    while i < (countChargesDIP2+countChargesDIP1):
        for j in range(len(ret1)):
            ret2.append(atom_xyz(ret1[j].symbol, ret1[j].coords))
        i += 1

    n = 0
    while n < countChargesPDIR:
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

def mergeChargeCoords(coords, charges, path):
    f = open(path, "w")
    for j in range(len(coords)):
        k = 0
        while k < 3:
            f.write(str(np.array(coords[j].coords[k])))
            f.write("   ")
            k += 1
        f.write(str(charges[j].charge))
        f.write("\n")
    f.close()

def mergeSymbolCoords(coords, symbols, path):
    f = open(path, "w")
    for j in range(len(coords)):
        k = 0
        f.write(symbols[j].symbol)
        while k < 3:
            f.write("   ")
            f.write(str(np.array(coords[j].coords[k])))
            k += 1
        f.write("\n")
    f.close()

'''Funktion output:
    Öffnet eine neue xyz-Datei und speichert dort
    die veränderten Koordinaten und Symbole ab.'''
def output1D(newcoord, fileName):
    f = open(fileName, "w")
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

# maybe cleanup later
def output2D(newcoord, fileName):
    f = open(fileName, "w")
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
def outputchargesDIP1D(newcoord, step):
    f = open("neuechargesDIP{}.txt".format(step), "w")
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

def outputchargesDIP2D(newcoord, step, step2):
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
def make_com(xyz1, xyz2, chargesDIP, chargesPDIR, char, method, basisSet, charMols):
    if basisSet == "1":
        basisSetString = "cc-pVDZ"
    if basisSet == "2":
        basisSetString = "STO-3G"
    if char == "Ja" and charMols == "Nein":
        if method == "2":
            head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
#p wB97XD/cc-pVDZ charge gfinput pop=full TD(nstates=8) 

Dimer

0 1
"""
        if method == "1" or method == "invalid":
            head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
# HF/{} charge

Dimer

0 1
""".format(basisSetString)
        if method == "3":
            head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
#p wB97XD/cc-pVDZ charge

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
        if method == "2":
            head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
#p wB97XD/cc-pVDZ gfinput pop=full TD(nstates=8) 

Dimer

0 1
"""
        if method == "1" or method == "invalid":
            head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
# HF/{}

Dimer

0 1
""".format(basisSetString)
        if method == "3":
            head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
#p wB97XD/cc-pVDZ

Dimer

0 1
"""
        if method == "4":
            head = """%NProcShared=8
%Mem=6GB
%chk=gs.chk
#p PM7

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
    Erstellt eine sh-Datei


    How to read the Job-Name:

    D1aP1CHS_08x_05y

    D1  1 DIP molecule
    a   DIP is rotated around a (b, c) axis
    P1  1 PDIR molecule (rotation currently not supported)
    C   calculated with Charges (or V for vacuum)
    H   Method (or w for wB97-XD or T for TD-cB97-XD or P for PM7)
    S   Basis set STO-3G (or c for cc-pVDZ)
    _   Space holder
    08  Step number
    x   axis of this stepping
    _   Space holder
   (05  optional second step number for 2D scans
    y   axis of this stepping)
    '''
def make_sh1D(nDIP, rotAxisDIP, nPDIR, chargeOrVacuum, method, basisSet, step1, axis1):
    head = """#PBS -l nodes=1:ppn={}
#PBS -l walltime=800:00:00
#PBS -l mem={}GB
#PBS -N D{}{}P{}{}{}{}_{}{}

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
    CV = ""
    if chargeOrVacuum == "charge":
        CV = "C"
    else:
        CV = "V"

    MET = ""
    PROC = "8"
    MEM = "8"
    if method == "1":
        MET = "H"
    if method == "2":
        MET = "T"
    if method == "3":
        MET = "w"
    if method == "4":
        MET = "P"
        PROC = "8"
        MEM = "8"

    BS = ""
    if basisSet == "0":
        BS = "-"
    if basisSet == "1":
        BS = "c"
    if basisSet == "2":
        BS = "S"

    with open("gauss.sh", "w") as out_file:
        out_file.write(head.format(PROC, MEM, nDIP, rotAxisDIP, nPDIR, CV, MET, BS, str(step1).zfill(2), axis1))

def make_sh2D(nDIP, rotAxisDIP, nPDIR, chargeOrVacuum, method, basisSet, step1, axis1, step2, axis2):
    head = """#PBS -l nodes=1:ppn={}
#PBS -l walltime=800:00:00
#PBS -l mem={}GB
#PBS -N D{}{}P{}{}{}{}_{}{}_{}{}

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
    CV = ""
    if chargeOrVacuum == "charge":
        CV = "C"
    else:
        CV = "V"

    MET = ""
    PROC = "8"
    MEM = "8"
    if method == "1":
        MET = "H"
    if method == "2":
        MET = "T"
    if method == "3":
        MET = "w"
    if method == "4":
        MET = "P"
        PROC = "4"
        MEM = "4"

    BS = ""
    if basisSet == "0":
        BS = "-"
    if basisSet == "1":
        BS = "c"
    if basisSet == "2":
        BS = "S"

    with open("gauss.sh","w") as out_file:
        out_file.write(head.format(PROC, MEM, nDIP, rotAxisDIP, nPDIR, CV, MET, BS, str(step1).zfill(2), axis1, str(step2).zfill(2), axis2))

def useroutput1D(char, DIP2, dup, dup2, verschiebung1, verschiebung2,
    xyz, shiftstart, shiftsize, shiftlength, wahlPDIR, axisPDIR, anglePDIR,
    wahlDIP, axisDIP, angleDIP, numberChargeLayersDIP, numberChargeLayersPDIR):
    f=open("Usereingaben.txt","w")
    f.write("Folgende Usereingaben wurden bei dieser Rechnung benutzt: ")
    f.write("\n\n")
    f.write("Wurden Ladungen berücksichtigt: %s" %char)
    f.write("\n")
    f.write("   DIP Charge Layer: %s" %numberChargeLayersDIP)
    f.write("\n")
    f.write("   PDIR Charge Layer: %s" %numberChargeLayersPDIR)
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
    f.write("Die DIP-Moleküle wurden um die Koordinate %s verschoben." %xyz)
    f.write("\n")
    f.write("Koordinate:")
    f.write("\n")
    f.write("   Startpunkt: %s" %shiftstart)
    f.write("\n")
    f.write("   Schrittlänge: %s" %shiftsize)
    f.write("\n")
    f.write("   Gesamtlänge: %s" %shiftlength)
    f.write("\n")
    f.write("Wurde der DIP-Kristall rotiert: %s" %wahlDIP)
    f.write("\n")
    f.write("   Achse: %s" %axisDIP)
    f.write("\n")
    f.write("   Winkel: %s" %angleDIP)
    f.write("\n")
    f.write("Wurde der PDIR-Kristall rotiert: %s" %wahlPDIR)
    f.write("\n")
    f.write("   Achse: %s" %axisPDIR)
    f.write("\n")
    f.write("   Winkel: %s" %anglePDIR)
    f.write("\n")
    f.close()

def useroutput2D(char, DIP2, dup, dup2, verschiebung1, verschiebung2,
    xyz1, xyz2, shiftstart, shiftstart2, shiftsize, shiftsize2, shiftlength,
    shiftlength2, wahl, axis, angle, abstandAuf3, numberChargeLayersDIP, numberChargeLayersPDIR):
    f = open("Usereingaben.txt","w")
    f.write("Folgende Usereingaben wurden bei dieser Rechnung benutzt: ")
    f.write("\n\n")
    f.write("Wurden Ladungen berücksichtigt: %s" %char)
    f.write("\n")
    f.write("   DIP Charge Layer: %s" %numberChargeLayersDIP)
    f.write("\n")
    f.write("   PDIR Charge Layer: %s" %numberChargeLayersPDIR)
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
    f.write("Abstand auf der dritten Koordinate: %s" %abstandAuf3)
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
            return "x"

def duplicateLayerDIP(original, offset):
    # Es ist egal ob man hier DIP1 oder DIP2 nimmt, denn es wird am Ende zusammengezählt.
    global countChargesDIP2
    newLayer = []
    for i in range(len(original)):
        countChargesDIP2 += 1
        currentCharge = copy.deepcopy(original[i])
        currentCharge.coords[2] += offset
        newLayer.append(charge_xyz(currentCharge.coords, currentCharge.charge))

    return newLayer



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

    axisDIP = ""
    angleDIP = ""

    axisPDIR = ""
    anglePDIR = ""

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

    #Winkel des DIP-Einheitskristalls
    alphaDIP = 90
    betaDIP = 92.416
    gammaDIP = 90

    #Differenz zwischen y- und b-Achse und z- und c-Achse
    difb = getbPDIR(gammaPDIR, lengthbPDIR)
    difc = getcPDIR(alphaPDIR, betaPDIR, gammaPDIR, lengthcPDIR)

    #Einlesen der Dateien
    file_input = readXYZ("output_first_DIP_molecule___with_opt_geo.xyz")
    file_input1 = readXYZ("output_second_DIP_molecule___with_opt_geo.xyz")
    # Choose which one you want:
    #file_input2 = readXYZ("PDIR_S0_1.xyz") #geometrieoptimierte Struktur
    file_input2 = readXYZ("PDIR_umsortiert.xyz") #lädt das Molekül in der Struktur des Kristalls
    file_input3 = readXYZ("PDIR_umsortiert.xyz") #lädt das Molekül in der Struktur des Kristalls

    #Initalisieren der Variablen mit Userinput
    mode = "invalid"
    while mode != "1" and mode != "2": # safety measure to ensure correct choice
        mode = input("Soll ein 1D oder 2D Scan ausgeführt werden? (1/2): ")

    sendCalculation = "invalid"
    while sendCalculation != "Ja" and sendCalculation != "Nein":
        sendCalculation = input("Soll die Rechnung am Ende abgeschickt werden? (Ja/Nein): ")

    calcMethod = "1"
    basisSet = "1" # 1: cc-pVDZ; 2: STO-3G;
    if sendCalculation =="Ja":
        calcMethod = "invalid"
        while calcMethod !="1" and calcMethod !="2" and calcMethod != "3" and calcMethod != "4":
            calcMethod = input("HF (1) oder TD-w-B97XD (2) oder w-B97XD (3) oder PM7 (4)? ")
        if calcMethod == "1":
            basisSet = input("Basissatz: cc-pVDZ (1) oder STO-3G (2)? ")
            while basisSet != "1" and basisSet != "2":
                basisSet = input("Basissatz: cc-pVDZ (1) oder STO-3G (2)? ")
    if calcMethod == "4":
        basisSet = "0"

    if mode == "1":
        DIP2 = input("Soll das 2. DIP-Molekül aus der Einheitszelle dazugeladen werden? (Ja, Nein): ")

        makeMoleculesFromCharges = "Nein"
        if calcMethod == "4":
            print("Bei PM7 werden Ladungen automatisch in Moleküle umgewandelt!")
            makeMoleculesFromCharges = "Ja"
        char = str(input("Sollen Ladungen berücksichtigt werden? (Ja, Nein): "))
        reduceChargeSize = "Nein"
        if char == "Ja":
            chargeOrVacuum = "charge"
        else:
            chargeOrVacuum = "vacuum"

        numberChargeLayersPDIR = 0
        numberChargeLayersDIP = 0
        if char == "Ja":
            numberChargeLayersPDIR = int(input("Wie viele Layer an PDIR Ladungen sollen berücksichtigt werden? "))
            numberChargeLayersDIP = int(input("Wie viele Layer an DIP Ladungen sollen berücksichtigt werden? "))
            if calcMethod != "4":
                makeMoleculesFromCharges = input("Sollen die Ladungen in Moleküle umgewandelt werden? (Ja, Nein): ")
            if makeMoleculesFromCharges == "Ja":
                reduceChargeSize = input("Soll die Ladungs-(Molekül-)Umgebung auf 3x3 reduziert werden? (Ja, Nein): ")

        dup = input("Soll DIP dupliziert werden (experimentell)? (Ja, Nein): ")
        dup2 = input("Soll PDIR dupliziert werden? (Ja, Nein): ")
        
        #Liest das 2. DIP-Molekül ein und hängt es an das erste an.
        nDIP = "1"
        if DIP2 == "Ja": # Use and char == "Nein" to center only on DIP1
            nDIP = "2"
            for i in range(len(file_input1)):
                file_input.append(atom_xyz(file_input1[i].symbol, file_input1[i].coords))



        '''Liest die Ladungen von DIP und PDIR ein, verschiebt sie,
            speichert sie in 2 Listen und macht centerofGeo von beiden'''
        CenterOfDIP = getCenterOfGeo(file_input)
        file_geo = moveToCenterofGeo(file_input)
        file_geo2 = moveToCenterofGeo(file_input2)
        file_geo3 = moveToCenterofGeo(file_input3) #PDIR charges crystal structure

        #Duplizieren sowohl von DIP als auch von PDIR
        verschiebung2 = "0"
        if dup == "Ja":
            verschiebung1 = input("In welche Richtung soll das zweite DIP-Molekül verschoben werden? (xp, xn, yp, yn, zp, zn): ")
            file_geo = duplicateDIP(file_geo, verschiebung1, lengthaDIP, lengthbDIP)
        
        nPDIR = "1"
        if dup2 == "Ja":
            verschiebung2 = input("In welche Richtung soll das zweite PDIR-Molekül verschoben werden? (ap, an, bp, bn, cp, cn): ") 
            file_geo2 = duplicatePDIR(file_geo2, lengthaPDIR, difb, difc, verschiebung2)
            nPDIR = "2"

        CenterOfPDIR = getCenterOfGeo(file_input2)

        # Vorgriff auf die Rotation, um zu wissen welche Charges ausgewählt werden müssen
        axisDIP = "0" # Wenn die 0 nicht überschrieben wird, wurde nicht gedreht.
        count =0
        wahlDIP = input("Willst du den DIP Kristall rotieren? (Ja, Nein): ")
        if wahlDIP == "Ja":
            axisDIP = input("Um welche Achse soll rotiert werden? (a, b, c): ")

        letteraxisDIP = axisDIP
        if wahlDIP == "Ja":
            if axisDIP == "a":
                axisDIP = [1,0,0]
                angleDIP = 180 - alphaDIP
                print("Vorsicht: y = z und z = y")
            elif axisDIP == "b":
                axisDIP = [0,1,0]
                angleDIP = 180 - betaDIP
                print("Vorsicht: x = z und z = x")
            elif axisDIP == "c":
                axisDIP = [0,0,1]
                angleDIP = 180 - gammaDIP
                print("Vorsicht: x = y und y = x")

        wahlPDIR = input("Willst du den PDIR-Kristall rotieren? (Ja, Nein): ")
        if wahlPDIR == "Ja":
            axisPDIR = input("Um welche Achse soll rotiert werden? (a, b, c): ")
            if axisPDIR == "a":
                axisPDIR = [1,0,0]
                anglePDIR = 180 - alphaPDIR

            elif axisPDIR == "b":
                axisPDIR = [0,1,0]
                anglePDIR = 180 - betaPDIR

            elif axisPDIR == "c":
                axisPDIR = [0,0,1]
                anglePDIR = 180 - gammaPDIR


        #Initalisieren der Variablen mit Userinput
        xyz = input("Gib die Koordinate an, die verschoben werden soll (x, y, z): ")
        
        shiftstart = float(input("Gib den Startpunkt der {}-Koordinate an: ".format(xyz)))
        shiftsize = float(input("Gib die Schrittweite der Verschiebung der {}-Koordinate an: ".format(xyz)))
        shiftlength = float(input("Gib den Endpunkt der Verschiebung ein: "))


        if char == "Ja":

            chargesDIP1 = readcharges("first_DIP_charges.txt")
            chargesDIP2 = readcharges("second_DIP_charges.txt")
            chargesPDIR = readcharges("PDIR_S0_1_charges.txt")
            mergeChargeCoords(file_input3, chargesPDIR, "Charges_Merged_PDIR.txt")
            chargesPDIR = readcharges("Charges_Merged_PDIR.txt")


            geo_chargesPDIR = moveToCenterofGeo(chargesPDIR)
            #geo_chargesDIP = moveToCenterofGeo(chargesDIP1)
            #geo_chargesDIP2 = moveToCenterofGeo(chargesDIP2)
            
        
            shiftchargesdip1 = getchargesDIP1(chargesDIP1, lengthaDIP, lengthbDIP, lengthcDIP, numberChargeLayersDIP, letteraxisDIP, xyz, reduceChargeSize) 
            shiftchargesdip2 = getchargesDIP2(chargesDIP2, lengthaDIP, lengthbDIP, lengthcDIP, numberChargeLayersDIP,
                DIP2, letteraxisDIP, xyz, reduceChargeSize) #DIP2: is the second DIP there?
            PDIRcharges = getchargesPDIR(geo_chargesPDIR, difb, lengthaPDIR, difc, letteraxisDIP, numberChargeLayersPDIR, verschiebung2, xyz, reduceChargeSize)
            #geo_chargesDIP = moveToCenterofGeo(chargesDIP1)
            #shiftchargesdip1 = moveToCenterofGeo(shiftchargesdip1)
            #shiftchargesdip2 = moveToCenterofGeo(shiftchargesdip2)
            #mergeSymbolCoords(chargesDIP2, file_input1, "NewDip2Molecule.txt")
            mergeChargeCoords(file_input1, chargesDIP2, "Charges_Merged_DIP2.txt")
            chargesDIP2 = readcharges("Charges_Merged_DIP2.txt")

            for i in range(len(shiftchargesdip2)):
                shiftchargesdip1.append(charge_xyz(shiftchargesdip2[i].coords, shiftchargesdip2[i].charge))
            # centerOfBothDIPs = getCenterOfGeo(shiftchargesdip1)
            if DIP2 == "Nein" and not(letteraxisDIP == "b" and (xyz == "z" or xyz == "y")):
                geo_chargesDIP = moveToCenterofGeo(shiftchargesdip1)
            
            if DIP2 == "Ja" or (letteraxisDIP == "b" and (xyz == "z" or xyz == "y")):
                for i in range(len(shiftchargesdip1)):
                    m = 0
                    while m < 3:
                        shiftchargesdip1[i].coords[m] -= CenterOfDIP[m]
                        m += 1
                geo_chargesDIP = copy.deepcopy(shiftchargesdip1)

            if dup2 == "Ja":
                for i in range(len(PDIRcharges)):
                    m = 0
                    while m < 3:
                        PDIRcharges[i].coords[m] -= CenterOfPDIR[m]
                        m += 1

            chargesBothDIP = copy.deepcopy(chargesDIP1)

            for i in range(len(chargesDIP2)):
                chargesBothDIP.append(charge_xyz(chargesDIP2[i].coords, chargesDIP2[i].charge))

            # mergeChargeCoords(file_input, chargesBothDIP, "MergedChargesCoordsBothDIP.txt")
            # chargesBothDIP = readcharges("MergedChargesCoordsBothDIP.txt")

            #centered_chargesBothDIP = moveToCenterofGeo(chargesBothDIP)
            if DIP2 == "Nein":
                geo_chargeDIPOriginal = moveToCenterofGeo(chargesDIP1) #The Charge that lies in the original DIP1 molecule
            if DIP2 == "Ja":
                for i in range(len(chargesDIP1)):
                    m = 0
                    while m < 3:
                        chargesBothDIP[i].coords[m] -= CenterOfDIP[m]
                        m += 1
                geo_chargeDIPOriginal=copy.deepcopy(chargesBothDIP)

            #geo_chargesDIP = shiftchargesdip1
            # duplicate the working first layer of DIP
            OneDirection = 1
            if letteraxisDIP == "a":
                OneDirection = -1

            OriginalOneDirection = OneDirection

            if numberChargeLayersDIP >= 2 or (letteraxisDIP == "a" and (xyz == "z" or xyz == "x")) or (letteraxisDIP == "b" and (xyz == "z" or xyz == "y")):
                DIPchargesSecondLayer = duplicateLayerDIP(geo_chargesDIP, OneDirection*lengthcDIP)
                if numberChargeLayersDIP >=3 or (letteraxisDIP == "a" and (xyz == "z" or xyz == "x")) or (letteraxisDIP == "b" and (xyz == "z" or xyz == "y")):
                    if (letteraxisDIP == "a" and (xyz == "z" or xyz == "x")) or (letteraxisDIP == "b" and (xyz == "z" or xyz == "y")):
                        OneDirection = -0.5*OneDirection # This modification makes the duplication function behave like expected when duplicating sideways
                    DIPchargesThirdLayer = duplicateLayerDIP(geo_chargesDIP, 2*OneDirection*lengthcDIP)
                    OneDirection = OriginalOneDirection
                    for i in range(len(DIPchargesThirdLayer)):
                        geo_chargesDIP.append(charge_xyz(DIPchargesThirdLayer[i].coords, DIPchargesThirdLayer[i].charge))

                for i in range(len(DIPchargesSecondLayer)):
                    geo_chargesDIP.append(charge_xyz(DIPchargesSecondLayer[i].coords, DIPchargesSecondLayer[i].charge))
                
                # Add in the charge directly over the original molecule
                Tempgeocharges = duplicateLayerDIP(geo_chargeDIPOriginal, OneDirection*lengthcDIP)
                for i in range(len(Tempgeocharges)):
                    geo_chargesDIP.append(charge_xyz(Tempgeocharges[i].coords, Tempgeocharges[i].charge))

                if numberChargeLayersDIP >= 3 or (letteraxisDIP == "a" and (xyz == "z" or xyz == "x")) or (letteraxisDIP == "b" and (xyz == "z" or xyz == "y")):
                    if (letteraxisDIP == "a" and (xyz == "z" or xyz == "x")) or (letteraxisDIP == "b" and (xyz == "z" or xyz == "y")):
                        OneDirection = -0.5*OneDirection
                    Tempgeocharges = duplicateLayerDIP(geo_chargeDIPOriginal, 2*OneDirection*lengthcDIP)
                    OneDirection = OriginalOneDirection
                    for i in range(len(Tempgeocharges)):
                        geo_chargesDIP.append(charge_xyz(Tempgeocharges[i].coords, Tempgeocharges[i].charge))
           
            geo_chargesPDIR = PDIRcharges

            

        #Erstellt Center of Geometry der Moleküle
        # CenterOfDIP = getCenterOfGeo(file_input)
        # file_geo = moveToCenterofGeo(file_input)  # Dip @ CoG
        #file_geo1 = moveToCenterofGeo(file_input1) #DIP2 @ CoG
        #This line has to be there to align the CoG of DIP1 and DIP2 to the CoG of the PDIR Dimer, no one knows why
        file_geo2 = moveToCenterofGeo(file_input2)
        

        ### Use this to center on DIP1 when using both DIPs
        # DIP2MoleculeFromCharges = readXYZ("NewDip2Molecule.txt")
        # # Move the DIP2 molecule from charges by the same amount as DIP1 was moved by moveToCenterOfGeo:
        # for i in range(len(DIP2MoleculeFromCharges)):
        #     m = 0
        #     while m < 3:
        #         DIP2MoleculeFromCharges[i].coords[m] -= CenterOfDIP[m]
        #         m += 1
        # # Add the corrected DIP2 to the DIP molecule
        # for oneAtom in DIP2MoleculeFromCharges:
        #     file_geo.append(atom_xyz(oneAtom.symbol, oneAtom.coords))
            
        if makeMoleculesFromCharges == "Ja":
            # Begin with PDIR:
            lenSinglePDIR = len(file_input2)
            j = 0
            for i in range(len(geo_chargesPDIR)):
                file_geo2.append(atom_xyz(file_input2[j].symbol, geo_chargesPDIR[i].coords))
                j += 1
                if j >= lenSinglePDIR:
                    j = 0
            # DIP:
            lenSingleDIP = len(file_input)
            j = 0
            for i in range(len(geo_chargesDIP)):
                file_geo.append(atom_xyz(file_input[j].symbol, geo_chargesDIP[i].coords))
                j += 1
                if j >= lenSingleDIP:
                    j = 0


        #Rotation des Kristalls
        if wahlDIP == "Ja":
            # From this point on axisDIP must be [ , , ]
            file_geo = product_matrix_vector(rot_axis_angle(axisDIP, angleDIP), file_geo)
            if char == "Ja":
                geo_chargesDIP = product_matrix_vector_charges(rot_axis_angle(axisDIP, angleDIP), geo_chargesDIP)
        

        #Rotation des Kristalls
        if wahlPDIR == "Ja":
            file_geo2 = product_matrix_vector(rot_axis_angle(axisPDIR, anglePDIR), file_geo2)
            if char == "Ja":
                geo_chargesPDIR = product_matrix_vector_charges(rot_axis_angle(axisPDIR, anglePDIR), geo_chargesPDIR)
        
        
        #Erstellt einen neuen Ordner, geht hinein und macht ihn zum Hauptordner
        new= input ("Ordnernamen angeben: ")
        os.mkdir(new)
        os.chdir(new)
        main_dir = os.getcwd()

        #Erstellt eine Datei, welche die Usereingaben beinhalten.

        useroutput1D(char, DIP2, dup, dup2, verschiebung1, verschiebung2, xyz, shiftstart,
            shiftsize, shiftlength, wahlPDIR, axisPDIR, anglePDIR, wahlDIP, axisDIP, angleDIP,
            numberChargeLayersDIP, numberChargeLayersPDIR)
        
        newdir="{}_{}_{}".format(str(step).zfill(2), xyz, chargeOrVacuum)
        os.mkdir(newdir)
        os.chdir(newdir)

        '''Berechnet den Startpunkt, speichert die Symbole und
            veränderten Koordinaten in eine Datei und liest sie wieder ein.'''
        shiftxyz = shiftcoord(xyz, file_geo, shiftstart)
        bothMolecules = copy.deepcopy(shiftxyz)
        for singleAtom in file_geo2:
            bothMolecules.append(atom_xyz(singleAtom.symbol, singleAtom.coords))

        output1D(bothMolecules, "neueKoordinatenBeide{}.xyz".format(step))
        output1D(shiftxyz, "neueKoordinaten{}.xyz".format(step))

        first_input = readXYZ("neueKoordinaten{}.xyz".format(step)) 

        '''Wenn Ladungen berücksichtigt werden:
            Verändert die Koordinaten der Ladungen zum Startpunkt, speichert sie in einer Datei
            liest sie ein und gibt eine xyz-Datei aus, in der die Ladungen und die beiden
            Moleküle DIP und PDIR vorhanden sind.'''
        if char == "Ja" and makeMoleculesFromCharges == "Nein":
            shiftchargexyz = shiftcoordcharges(xyz, geo_chargesDIP, shiftstart)
            outputchargesDIP1D(shiftchargexyz, step)

            first_charges = readcharges("neuechargesDIP{}.txt".format(step))

            DIPPDIR = getonefile(file_geo2, first_input)
            # first_input: Um variable verschobenes, vorher CoG aligntes DIP (echtes mol)
            # filegeo2: CoG aligned PDIR (echtes mol)
            # DIPPDIR: first_input + file_geo2 (echtes mol)
            # First_charges: Um var. verschobene, vorher CoG alignte DIP-Charges
            # geo_chargesPDIR: CoG aligned PDIR-Charges, nicht verschoben
            outputeverythingDIPPDIR1D(first_input, file_geo2, DIPPDIR, first_charges,
                geo_chargesPDIR, step)

        #Erstellen der com- und sh-Datei
        make_com(first_input, file_geo2, first_charges, geo_chargesPDIR, char, calcMethod, basisSet, makeMoleculesFromCharges)
        make_sh1D(nDIP, letteraxisDIP, nPDIR, chargeOrVacuum, calcMethod, basisSet, step, xyz)

        #Abschicken der Rechnung
        if sendCalculation == "Ja":
            subprocess.call('qsub -V gauss.sh', shell=True)

        #Schleife, die die ebengenannten Schritte wiederholt.
        i = shiftstart
        while i < shiftlength:

            new_input = readXYZ("neueKoordinaten%s.xyz" %step)

            if char == "Ja" and makeMoleculesFromCharges == "Nein":
                    new_charges = readcharges("neuechargesDIP{}.txt".format(step))

            step = step + 1

            os.chdir(main_dir)

            newdir="{}_{}_{}".format(str(step).zfill(2), xyz, chargeOrVacuum)
            os.mkdir(newdir)
            os.chdir(newdir)

            shiftxyzstep = shiftcoord(xyz, new_input, shiftsize)
            bothMoleculesstep = copy.deepcopy(shiftxyzstep)
            for singleAtom in file_input2:
                bothMoleculesstep.append(atom_xyz(singleAtom.symbol, singleAtom.coords))

            output1D(bothMoleculesstep, "neueKoordinatenBeide{}.xyz".format(step))
            output1D(shiftxyzstep, "neueKoordinaten{}.xyz".format(step))
            if char == "Ja" and makeMoleculesFromCharges == "Nein":
                shiftchargexyz = shiftcoordcharges(xyz, new_charges, shiftsize)
                outputchargesDIP1D(shiftchargexyz, step)
                new_charges = readcharges("neuechargesDIP{}.txt".format(step))
                DIPPDIR = getonefile(file_geo2, new_input)
                outputeverythingDIPPDIR1D(new_input, file_geo2, DIPPDIR, new_charges,
                    geo_chargesPDIR, step)
                
            

            make_com(new_input, file_geo2, new_charges, geo_chargesPDIR, char, calcMethod, basisSet, makeMoleculesFromCharges)
            make_sh1D(nDIP, letteraxisDIP, nPDIR, chargeOrVacuum, calcMethod, basisSet, step, xyz)
            
            if sendCalculation == "Ja":                
                subprocess.call('qsub -V gauss.sh', shell=True)
                                    
            i += shiftsize


    '''Begin 2D Mode




    Visual aid for better visibility.






    '''

    if mode == "2":
        DIP2 = input("Soll das 2. DIP-Molekül aus der Einheitszelle dazugeladen werden? (Ja, Nein): ")
        char = "Nein"
        makeMoleculesFromCharges = "Nein"
        if calcMethod == "4":
            print("Bei PM7 werden Ladungen automatisch in Moleküle umgewandelt!")
            makeMoleculesFromCharges = "Ja"
        char = str(input("Sollen Ladungen berücksichtigt werden? (Ja, Nein): "))
        reduceChargeSize = "Nein"
        if char == "Ja":
            chargeOrVacuum = "charge"
        else:
            chargeOrVacuum = "vacuum"

        numberChargeLayersPDIR = 0
        numberChargeLayersDIP = 0
        if char == "Ja":
            numberChargeLayersPDIR = int(input("Wie viele Layer an PDIR Ladungen sollen berücksichtigt werden? "))
            numberChargeLayersDIP = int(input("Wie viele Layer an DIP Ladungen sollen berücksichtigt werden? "))
            if calcMethod != "4":
                makeMoleculesFromCharges = input("Sollen die Ladungen in Moleküle umgewandelt werden? (Ja, Nein): ")
            if makeMoleculesFromCharges == "Ja":
                reduceChargeSize = input("Soll die Ladungs-(Molekül-)Umgebung auf 3x3 reduziert werden? (Ja, Nein): ")
        dup = input("Soll DIP dupliziert werden (experimentell)? (Ja, Nein): ")
        dup2 = input("Soll PDIR dupliziert werden? (Ja, Nein): ")
        
        #Liest das 2. DIP-Molekül ein und hängt es an das erste an.
        nDIP = "1"
        if DIP2 == "Ja":
            nDIP = "2"
            for i in range(len(file_input1)):
                file_input.append(atom_xyz(file_input1[i].symbol, file_input1[i].coords))

        #Erstellt Center of Geometry der Moleküle
        CenterOfDIP = getCenterOfGeo(file_input)
        file_geo = moveToCenterofGeo(file_input)  # Dip @ CoG
        #file_geo1 = moveToCenterofGeo(file_input1) #DIP2 @ CoG
        file_geo2 = moveToCenterofGeo(file_input2)
        file_geo3 = moveToCenterofGeo(file_input3) #PDIR charges crystal structure

        
        #Duplizieren sowohl von DIP als auch von PDIR
        verschiebung2 = "0"
        if dup == "Ja":
            verschiebung1 = input("In welche Richtung soll das zweite DIP-Molekül verschoben werden? (xp, xn, yp, yn, zp, zn):")
            file_geo = duplicateDIP(file_geo, verschiebung1, lengthaDIP, lengthbDIP)
            # for i in range(len(file_geo_verschoben)): # No longer needed, probably erroneous
            #     file_geo.append(atom_xyz(file_geo_verschoben[i].symbol, file_geo_verschoben[i].coords))

        nPDIR = "1"
        if dup2 == "Ja":
            verschiebung2 = input("In welche Richtung soll das zweite PDIR-Molekül verschoben werden? (ap, an, bp, bn, cp, cn):") 
            file_geo2 = duplicatePDIR(file_geo2, lengthaPDIR, difb, difc, verschiebung2)
            nPDIR = "2"
            # for i in range(len(file_geo2_verschoben)): # No longer needed, probably erroneous
            #     file_geo2.append(atom_xyz(file_geo2_verschoben[i].symbol, file_geo2_verschoben[i].coords))

        #file_geo3 = moveToCenterofGeo(file_input3) #PDIR charges crystal structure
        CenterOfPDIR = getCenterOfGeo(file_input2)

        # Vorgriff auf die Rotation, um zu wissen welche Charges ausgewählt werden müssen
        axisDIP = "0" # Wenn die 0 nicht überschrieben wird, wurde nicht gedreht.
        count =0
        wahlDIP = input("Willst du den DIP Kristall rotieren? (Ja, Nein): ")
        if wahlDIP == "Ja":
            axisDIP = input("Um welche Achse soll rotiert werden? (a, b, c): ")

        letteraxisDIP = axisDIP
        if wahlDIP == "Ja":
            if axisDIP == "a":
                axisDIP = [1,0,0]
                angleDIP = 180 - alphaDIP
                print("Vorsicht: y = z und z = y")
            elif axisDIP == "b":
                axisDIP = [0,1,0]
                angleDIP = 180 - betaDIP
                print("Vorsicht: x = z und z = x")
            elif axisDIP == "c":
                axisDIP = [0,0,1]
                angleDIP = 180 - gammaDIP
                print("Vorsicht: x = y und y = x")

        wahlPDIR = input("Willst du den PDIR-Kristall rotieren? (Ja, Nein): ")
        if wahlPDIR == "Ja":
            axisPDIR = input("Um welche Achse soll rotiert werden? (a, b, c): ")
            if axisPDIR == "a":
                axisPDIR = [1,0,0]
                anglePDIR = 180 - alphaPDIR

            elif axisPDIR == "b":
                axisPDIR = [0,1,0]
                anglePDIR = 180 - betaPDIR

            elif axisPDIR == "c":
                axisPDIR = [0,0,1]
                anglePDIR = 180 - gammaPDIR

        #Initalisieren weiterer Varibalen
        xyz1 = input("Gib die 1. Koordinate an, die verschoben werden soll (x, y, z): ")
        xyz2 = input("Gib die 2. Koordinate an, die verschoben werden soll (x, y, z): ")
        distance = float(input("Gib den Abstand auf der {}-Achse an: ".format(thirdAxis(xyz1, xyz2))))

        shiftstart = float(input("Gib den Startpunkt der ersten Koordinate an: "))
        shiftsize = float(input("Gib die Schrittweite der Verschiebung der ersten Koordinate an: "))
        shiftlength = float(input("Gib den Endpunkt der erste Koordinate an: "))
        shiftstart2 = float(input("Gib den Startpunkt der zweiten Koordinate an: "))
        shiftsize2 = float(input("Gib die Schrittweite der Verschiebung der zweiten Koordinate an: "))
        shiftlength2 = float(input("Gib den Endpunkt der zweiten Koordinate an: "))
        
        if char == "Ja":

            chargesDIP1 = readcharges("first_DIP_charges.txt")
            chargesDIP2 = readcharges("second_DIP_charges.txt")
            chargesPDIR = readcharges("PDIR_S0_1_charges.txt")
            mergeChargeCoords(file_input3, chargesPDIR, "Charges_Merged_PDIR.txt")
            chargesPDIR = readcharges("Charges_Merged_PDIR.txt")


            geo_chargesPDIR = moveToCenterofGeo(chargesPDIR)
            #geo_chargesDIP = moveToCenterofGeo(chargesDIP1)
            #geo_chargesDIP2 = moveToCenterofGeo(chargesDIP2)
            
        
            shiftchargesdip1 = getchargesDIP1(chargesDIP1, lengthaDIP, lengthbDIP, lengthcDIP, numberChargeLayersDIP, letteraxisDIP, thirdAxis(xyz1, xyz2), reduceChargeSize) 
            shiftchargesdip2 = getchargesDIP2(chargesDIP2, lengthaDIP, lengthbDIP, lengthcDIP, numberChargeLayersDIP,
                DIP2, letteraxisDIP, thirdAxis(xyz1, xyz2), reduceChargeSize) #DIP2: is the second DIP there?
            PDIRcharges = getchargesPDIR(geo_chargesPDIR, difb, lengthaPDIR, difc, letteraxisDIP, numberChargeLayersPDIR, verschiebung2, thirdAxis(xyz1, xyz2), reduceChargeSize)
            #geo_chargesDIP = moveToCenterofGeo(chargesDIP1)
            #shiftchargesdip1 = moveToCenterofGeo(shiftchargesdip1)
            #shiftchargesdip2 = moveToCenterofGeo(shiftchargesdip2)
            #mergeSymbolCoords(chargesDIP2, file_input1, "NewDip2Molecule.txt")
            mergeChargeCoords(file_input1, chargesDIP2, "Charges_Merged_DIP2.txt")
            chargesDIP2 = readcharges("Charges_Merged_DIP2.txt")

            xyz3 = thirdAxis(xyz1, xyz2)

            for i in range(len(shiftchargesdip2)):
                shiftchargesdip1.append(charge_xyz(shiftchargesdip2[i].coords, shiftchargesdip2[i].charge))
            # centerOfBothDIPs = getCenterOfGeo(shiftchargesdip1)
            if DIP2 == "Nein" and not(letteraxisDIP == "b" and (xyz3 == "z" or xyz3 == "y")):
                geo_chargesDIP = moveToCenterofGeo(shiftchargesdip1)
            
            if DIP2 == "Ja" or (letteraxisDIP == "b" and (xyz3 == "z" or xyz3 == "y")):
                for i in range(len(shiftchargesdip1)):
                    m = 0
                    while m < 3:
                        shiftchargesdip1[i].coords[m] -= CenterOfDIP[m]
                        m += 1
                geo_chargesDIP = copy.deepcopy(shiftchargesdip1)

            if dup2 == "Ja":
                for i in range(len(PDIRcharges)):
                    m = 0
                    while m < 3:
                        PDIRcharges[i].coords[m] -= CenterOfPDIR[m]
                        m += 1

            chargesBothDIP = copy.deepcopy(chargesDIP1)

            for i in range(len(chargesDIP2)):
                chargesBothDIP.append(charge_xyz(chargesDIP2[i].coords, chargesDIP2[i].charge))

            # mergeChargeCoords(file_input, chargesBothDIP, "MergedChargesCoordsBothDIP.txt")
            # chargesBothDIP = readcharges("MergedChargesCoordsBothDIP.txt")

            #centered_chargesBothDIP = moveToCenterofGeo(chargesBothDIP)
            if DIP2 == "Nein":
                geo_chargeDIPOriginal = moveToCenterofGeo(chargesDIP1) #The Charge that lies in the original DIP1 molecule
            if DIP2 == "Ja":
                for i in range(len(chargesDIP1)):
                    m = 0
                    while m < 3:
                        chargesBothDIP[i].coords[m] -= CenterOfDIP[m]
                        m += 1
                geo_chargeDIPOriginal=copy.deepcopy(chargesBothDIP)

            #geo_chargesDIP = shiftchargesdip1
            # duplicate the working first layer of DIP
            OneDirection = 1
            if letteraxisDIP == "a":
                OneDirection = -1
            OriginalOneDirection = OneDirection
            

            if numberChargeLayersDIP >= 2 or (letteraxisDIP == "a" and (xyz3 == "z" or xyz3 == "x")) or (letteraxisDIP == "b" and (xyz3 == "z" or xyz3 == "y")):
                DIPchargesSecondLayer = duplicateLayerDIP(geo_chargesDIP, OneDirection*lengthcDIP)
                if numberChargeLayersDIP >=3 or (letteraxisDIP == "a" and (xyz3 == "z" or xyz3 == "x")) or (letteraxisDIP == "b" and (xyz3 == "z" or xyz3 == "y")):
                    if letteraxisDIP == "a" and (xyz3 == "z" or xyz3 == "x") or (letteraxisDIP == "b" and (xyz3 == "z" or xyz3 == "y")):
                        OneDirection = -0.5*OneDirection
                    DIPchargesThirdLayer = duplicateLayerDIP(geo_chargesDIP, 2*OneDirection*lengthcDIP)
                    OneDirection = OriginalOneDirection
                    for i in range(len(DIPchargesThirdLayer)):
                        geo_chargesDIP.append(charge_xyz(DIPchargesThirdLayer[i].coords, DIPchargesThirdLayer[i].charge))

                for i in range(len(DIPchargesSecondLayer)):
                    geo_chargesDIP.append(charge_xyz(DIPchargesSecondLayer[i].coords, DIPchargesSecondLayer[i].charge))
                
                # Add in the charge directly over the original molecule
                Tempgeocharges = duplicateLayerDIP(geo_chargeDIPOriginal, OneDirection*lengthcDIP)
                for i in range(len(Tempgeocharges)):
                    geo_chargesDIP.append(charge_xyz(Tempgeocharges[i].coords, Tempgeocharges[i].charge))

                if numberChargeLayersDIP >= 3 or (letteraxisDIP == "a" and (xyz3 == "z" or  xyz3 == "x")) or (letteraxisDIP == "b" and (xyz3 == "z" or xyz3 == "y")):
                    if letteraxisDIP == "a" and (xyz3 == "z" or xyz3 == "x") or (letteraxisDIP == "b" and (xyz3 == "z" or xyz3 == "y")):
                        OneDirection = -0.5*OneDirection
                    Tempgeocharges = duplicateLayerDIP(geo_chargeDIPOriginal, 2*OneDirection*lengthcDIP)
                    OneDirection = OriginalOneDirection
                    for i in range(len(Tempgeocharges)):
                        geo_chargesDIP.append(charge_xyz(Tempgeocharges[i].coords, Tempgeocharges[i].charge))
           
            geo_chargesPDIR = PDIRcharges

            
        #This line has to be there to align the CoG of DIP1 and DIP2 to the CoG of the PDIR Dimer, no one knows why
        file_geo2 = moveToCenterofGeo(file_input2)

        if makeMoleculesFromCharges == "Ja":
            # Begin with PDIR:
            lenSinglePDIR = len(file_input2)
            j = 0
            for i in range(len(geo_chargesPDIR)):
                file_geo2.append(atom_xyz(file_input2[j].symbol, geo_chargesPDIR[i].coords))
                j += 1
                if j >= lenSinglePDIR:
                    j = 0
            # DIP:
            lenSingleDIP = len(file_input)
            j = 0
            for i in range(len(geo_chargesDIP)):
                file_geo.append(atom_xyz(file_input[j].symbol, geo_chargesDIP[i].coords))
                j += 1
                if j >= lenSingleDIP:
                    j = 0

        #Rotation des Kristalls
        if wahlDIP == "Ja":
            file_geo = product_matrix_vector(rot_axis_angle(axisDIP, angleDIP), file_geo)
            if char == "Ja":
                geo_chargesDIP = product_matrix_vector_charges(rot_axis_angle(axisDIP, angleDIP), geo_chargesDIP)

        #Rotation des Kristalls
        if wahlPDIR == "Ja":
            file_geo2 = product_matrix_vector(rot_axis_angle(axisPDIR, anglePDIR), file_geo2)
            if char == "Ja":
                geo_chargesPDIR = product_matrix_vector_charges(rot_axis_angle(axisPDIR, anglePDIR), geo_chargesPDIR)


        #Erstellt einen neuen Ordner, geht hinein und macht ihn zum Hauptordner
        new= input("Ordnername: ")
        os.mkdir(new)
        os.chdir(new)
        main_dir = os.getcwd()

        #Erstellt eine Datei, welche die Usereingaben beinhalten.
        useroutput2D(char, DIP2, dup, dup2, verschiebung1, verschiebung2,
            xyz1, xyz2, shiftstart, shiftstart2, shiftsize, shiftsize2,
            shiftlength, shiftlength2, wahlDIP, axisDIP, angleDIP, distance,
            numberChargeLayersDIP, numberChargeLayersPDIR)

        #Erstellt einen neuen Ordner und geht hinein
        newdir="{}_{}_{}".format(str(step).zfill(2), xyz1, chargeOrVacuum)
        os.mkdir(newdir)
        os.chdir(newdir)

        #Erstellt einen neuen Ordner und geht hinein
        newdir2 = "{}_{}_{}".format(str(step2).zfill(2), xyz2, chargeOrVacuum)
        os.mkdir(newdir2)
        os.chdir(newdir2)

        '''Berechnet den Startpunkt, speichert die Symbole und
            veränderten Koordinaten in eine Datei und liest sie wieder ein.'''
        shiftxyz1 = shiftcoord(xyz1, file_geo, shiftstart)
        shiftz = shiftcoord(thirdAxis(xyz1, xyz2), shiftxyz1, distance)
        shiftall = shiftcoord(xyz2, shiftz, shiftstart2)

        shiftallboth = copy.deepcopy(shiftall)
        for singleAtom in file_geo2:
                shiftallboth.append(atom_xyz(singleAtom.symbol, singleAtom.coords))

        output2D(shiftallboth, "neueKoordinatenBeide{}_{}.xyz".format(step, step2))
        output2D(shiftall, "neueKoordinaten{}_{}.xyz".format(step, step2))

        first_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))

        '''Wenn Ladungen berücksichtigt werden:
            Verändert die Koordinaten der Ladungen zum Startpunkt, speichert sie in einer Datei
            liest sie ein und gibt eine xyz-Datei aus, in der die Ladungen und die beiden
            Moleküle DIP und PDIR vorhanden sind.'''
        if char == "Ja" and makeMoleculesFromCharges == "Nein":
            shiftchargexyz1 = shiftcoordcharges(xyz1, geo_chargesDIP, shiftstart)
            shiftchargez = shiftcoordcharges(thirdAxis(xyz1, xyz2), shiftchargexyz1, distance)
            shiftallcharges = shiftcoordcharges(xyz2, shiftchargez, shiftstart2)

            outputchargesDIP2D(shiftallcharges, step, step2)

            first_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
            
            DIPPDIR = getonefile(file_geo2, first_input)

            outputeverythingDIPPDIR2D(first_input, file_geo2, DIPPDIR, first_charges, geo_chargesPDIR, step, step2)

        #Erstellen der com- und sh-Datei
        make_com(first_input, file_geo2, first_charges, geo_chargesPDIR, char, calcMethod, basisSet, makeMoleculesFromCharges)
        make_sh2D(nDIP, letteraxisDIP, nPDIR, chargeOrVacuum, calcMethod, basisSet, step, xyz1, step2, xyz2)

        #Abschicken der Rechnung
        if sendCalculation == "Ja":
            subprocess.call('qsub -V gauss.sh', shell = True)

        #Schleife, die die ebengenannten Schritte wiederholt.
        i = shiftstart
        while i <= shiftlength:

            j = shiftstart2
            while j < shiftlength2:

                new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))

                if char == "Ja" and makeMoleculesFromCharges == "Nein":
                    new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
                    
                step2 = step2 + 1

                os.chdir(main_dir)
                os.chdir(newdir)

                newdir3 = "{}_{}_{}".format(str(step2).zfill(2), xyz2, chargeOrVacuum)
                os.mkdir(newdir3)
                os.chdir(newdir3)

                shifty = shiftcoord(xyz2, new_input, shiftsize2)

                shiftallbothstepy = copy.deepcopy(shifty)
                for singleAtom in file_input2:
                    shiftallbothstepy.append(atom_xyz(singleAtom.symbol, singleAtom.coords))

                output2D(shiftallbothstepy, "neueKoordinatenBeide{}_{}.xyz".format(step, step2))

                output2D(shifty, "neueKoordinaten{}_{}.xyz".format(step, step2))

                if char == "Ja" and makeMoleculesFromCharges == "Nein":
                    shiftchargey = shiftcoordcharges(xyz2, new_charges, shiftsize2)

                    outputchargesDIP2D(shiftchargey, step, step2)

                    new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))
                    new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
                    DIPPDIR = getonefile(file_geo2, new_input)

                    outputeverythingDIPPDIR2D(new_input, file_geo2, DIPPDIR, new_charges, geo_chargesPDIR, step, step2)
                
     
                make_com(new_input, file_geo2, new_charges, geo_chargesPDIR, char, calcMethod, basisSet, makeMoleculesFromCharges)
                make_sh2D(nDIP, letteraxisDIP, nPDIR, chargeOrVacuum, calcMethod, basisSet, step, xyz1, step2, xyz2)

                if sendCalculation == "Ja":
                    subprocess.call('qsub -V gauss.sh', shell = True)

                j = j + shiftsize2
            
            if i < shiftlength:
                step2 = 0

                os.chdir(main_dir)
                os.chdir(newdir)
                os.chdir(newdir2)

                new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))

                if char == "Ja" and makeMoleculesFromCharges == "Nein":
                    new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))
                    new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))

                step = step + 1
            
                os.chdir(main_dir)

                newdir = "{}_{}_{}".format(str(step).zfill(2), xyz1, chargeOrVacuum)
                os.mkdir(newdir)
                os.chdir(newdir)

                newdir5 = "{}_{}_{}".format(str(step2).zfill(2), xyz2, chargeOrVacuum)
                os.mkdir(newdir5)
                os.chdir(newdir5)

                shiftx = shiftcoord(xyz1, new_input, shiftsize)

                shiftallbothstepx = copy.deepcopy(shiftx)
                for singleAtom in file_input2:
                    shiftallbothstepx.append(atom_xyz(singleAtom.symbol, singleAtom.coords))

                output2D(shiftallbothstepx, "neueKoordinatenBeide{}_{}.xyz".format(step, step2))

                output2D(shiftx, "neueKoordinaten{}_{}.xyz".format(step, step2))

                if char == "Ja" and makeMoleculesFromCharges == "Nein":
                    shiftchargex = shiftcoordcharges(xyz1, new_charges, shiftsize)

                    outputchargesDIP2D(shiftchargex, step, step2)

                    new_input = readXYZ("neueKoordinaten{}_{}.xyz".format(step, step2))
                    new_charges = readcharges("neuechargesDIP{}_{}.txt".format(step, step2))
                    DIPPDIR = getonefile(file_geo2, new_input)
                    
                    outputeverythingDIPPDIR2D(new_input, file_geo2, DIPPDIR, new_charges, geo_chargesPDIR, step, step2)

                make_com(new_input, file_geo2, new_charges, geo_chargesPDIR, char, calcMethod, basisSet, makeMoleculesFromCharges)
                make_sh2D(nDIP, letteraxisDIP, nPDIR, chargeOrVacuum, calcMethod, basisSet, step, xyz1, step2, xyz2)

                if sendCalculation == "Ja":
                    subprocess.call('qsub -V gauss.sh', shell = True)

            else:
                break
       
            i = i+shiftsize


    
if __name__=="__main__": main()   





