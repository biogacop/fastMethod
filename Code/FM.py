#!/usr/bin/env python
# coding: utf-8

import cobra
import time
import numpy as np
import more_itertools as mi
import multiprocessing as mp
import random
import gc
import pickle
from cobra import Model, Reaction, Metabolite

def identificaTipo(nombre):
    if "(" in nombre:
        indice1=nombre.index("(")
        nombre2=nombre[indice1+1:]
        indice2=nombre2.index(")")
        nombre3=nombre2[:indice2] 
        if " and " in nombre3:
            tipo="SOP"
        if " or " in nombre3:
            tipo="POS"
    else:
        tipo="SOP"
    return(tipo)

def toPOS(nombre,tope=6):
    productos=set()    
    tipo=identificaTipo(nombre)
    if tipo=="SOP":
        sumandos=nombre.split(" or ")
        for sumando in sumandos:
            if " and " in sumando:
                casos=sumando.split(" and ")
                miCaso=set()
                for caso in casos:
                    caso2=caso.replace("(","")
                    caso3=caso2.replace(")","")
                    caso4=caso3.replace(" ","")
                    miCaso.add(caso4)
                productos.add(frozenset(miCaso))
            else:
                productos.add(frozenset([sumando]))
        CS=set([frozenset()])
        for producto in productos:
            CS2=set()
            for C in CS:
                if len(set(C).intersection(producto))==0:
                    for a in producto:
                        C2=set(C.copy())
                        C2.add(a)
                        if len(C2)<=tope:
                            CS2.add(frozenset(C2))
                else:
                    CS2.add(C)
            CS=CS2
        CS3=set()
        for C in CS:
            bien=True
            for C2 in CS:
                if not C==C2 and C2.issubset(C):
                    bien=False
                    break
            if bien:
                CS3.add(C)
    else:
        CS2=set()
        productos = nombre.split(" and ")
        for producto in productos:
            CS=set()
            componentes=producto.split(" or ")
            for C in componentes:
                C2=C.replace("(","")
                C3=C2.replace(")","")
                C4=C3.replace(" ","")
                CS.add(C4)
            CS2.add(frozenset(CS))
        CS3=set()
        for C in CS2:
            bien=True
            for C2 in CS2:
                if not C==C2 and C2.issubset(C):
                    bien=False
                    break
            if bien:
                CS3.add(C)
    return(CS3)

def asociar(modelo,tope=6,genes=set()):
    total=0
    asociados=dict()
    reacciones=dict()
    eliminar=set()

    for k in range(tope+1):
        for num in range(len(modelo.reactions)):
            asociados[(modelo.reactions[num],k)]=set()        
    for num in range(len(modelo.reactions)):
        nombre=modelo.reactions[num].gene_reaction_rule
        CS3=toPOS(nombre,tope)
        k=0
        while k<=tope:
            mas=0
            nivel=set([C for C in CS3 if len(C)==k and not C==frozenset({''})])
            if len(nivel)>0:
                mios=set()
                for caso in nivel:
                    mio=set()
                    for a in caso:
                        mio.add(a)
                    l=len(mio)
                    mio=mio.difference(genes)
                    if l>0 and len(mio)==0:
                        eliminar.add(modelo.reactions[num].id)
                    asociados[(modelo.reactions[num],len(mio))].add(frozenset(mio))

                    mios.add(frozenset(mio))

                for caso in mios:
                    if not caso in reacciones:
                        reacciones[caso]=set([modelo.reactions[num]])
                    else:
                        reacciones[caso].add(modelo.reactions[num])
                mas+=len(nivel)
            total+=mas
            k+=1    
        
    return(asociados,reacciones,eliminar)
    
    
def soporteG(soporte,tope):
    sopG1=set(frozenset(case) for j in soporte for k in range(1,tope+1) for case in asociados[(modelo.reactions[j],k)])
    miTope=max(set(len(case) for case in sopG1))
    menores=[set()]
    for l in range(miTope):
        menores.append(set(case for case in sopG1 if len(case)==l))
    sopG2=set()
    for caso in sopG1:
        if caso in powers:
            po=powers[caso]
        else:
            po=set([frozenset(a) for a in mi.powerset(caso) if len(a)>0])
            powers[caso]=po
        misMenores=set()
        for l in range(len(caso)):
            misMenores.update(menores[l])
        if po.isdisjoint(misMenores):
            sopG2.add(caso)
    return(sopG2)

def chequea(inicio,fin):
    respuestas=set()
    for i in range(inicio,fin):
        if i<len(finales):
            C=finales[i]
            vale=True
            if not C.isdisjoint(CSs1):
                vale = False
            else:
                if C in powers:
                    po=powers[C]
                else:
                    po=set([frozenset(a) for a in mi.powerset(C) if len(a)>0])
                    powers[C]=po
                if not po.isdisjoint(todosMCS):
                    vale=False
                if vale:
                    for sop in listaSoportes[indice[C]:]:                    
                        if po.isdisjoint(sopG[sop]):           
                            vale=False
                            break
            if vale:
                respuestas.add(C)
    return(respuestas)


def potencia(caso):
    if caso in powers:
        return(powers[caso])
    else:
        return(set([frozenset(a) for a in mi.powerset(caso) if len(a)>0]))
    
def reduceLista2():
    tipos=[]
    for i in range(tope+1):
        tipos.append(set())
    for i in range(tope+1):
        tipos[i]=set(C for C in lista2 if len(C)==i)
    validos=tipos[1]
    subs=tipos[1]
    for i in range(2,tope+1):
        validos.update(set(C for C in tipos[i] if potencia(C).isdisjoint(subs)))
        subs.update(validos)
    return(validos)

def CSV2(_modelo,tope):
    global todosMCS, misSoportes, sopG, powers, modelo, miLista2, listaL2, finales
    global CSs1, nuevos, listaSoportes, indice,iMio, lista2, listaGenes
    
    modelo=_modelo
    listaMCS=[]
    listaMCS.append(set())    
    listaCandidatos=[]
    listaCandidatos.append(set([frozenset()]))
    todosMCS=set()
    cutsets=[set()]
    candidatos=[set()]

    
    for i in range(1,tope+1):
        if i==1:      
            with modelo:
                sol=modelo.optimize()
                sop=frozenset([j for j in range(len(modelo.reactions)) if abs(sol.fluxes.iloc[j])>10**-12])
                sopG2=soporteG(sop,1)

            miL2=set()
            listaGenes=[g for C in sopG2 for g in C]
            misSoportes.add(sop)
            
            k=0
            argus=[]
            num=miNum0
            cpus=cpus0
            while k*num<=len(listaGenes):
                argus.append((k*num,(k+1)*num+1))
                k+=1

            with mp.Pool(processes=cpus) as mp_pool:
                results=mp_pool.starmap(analizaGenes,argus)

            CS0=set()
            for result in results:
                CS0.update(result[0])
                misSoportes.update(result[1])
                
            todosMCS.update(CS0)    
            cutsets.append(CS0)
            
            for sop in misSoportes:
                sopG[sop]=soporteG(sop,tope)
            

            CSs1=set([g for caso in CS0 for g in caso])
            
       
        ############################## Final Fase 1 ############################################
    
        else:
            iMio=i
            
            CSs=set()
            for C in cutsets:
                CSs.update(C)

            for sop in misSoportes:
                sopG[sop].difference_update(CSs)
            ############################## Parte 1.- Berge #################################### 
            CST2=set()
            lista=set([frozenset()])
            descartados=set()
            cont2=0
            long=dict()
            indice=dict()
            analizados=set()
            
            for sop in misSoportes:
                long[sop]=len(sopG[sop])
            sortedSupports= sorted(long.items(), key=lambda x:x[1])
            finales=set()
            listaSoportes=[sortedSupports[t][0] for t in range(len(misSoportes))]
            
            for t in range(len(misSoportes)):
                ################ Para cada soporte #################
                item=sortedSupports[t]
                sop=item[0]
                
                lista2=set()
                for caso in lista:
                    if caso in powers:
                        po=powers[caso]
                    else:
                        po=set([frozenset(a) for a in mi.powerset(caso) if len(a)>0])
                        powers[caso]=po

                    if len(po.intersection(sopG[sop]))>0:
                        lista2.add(caso)
                    else:
                        for caso2 in sopG[sop]:
                            if len(caso.union(caso2))<iMio:
                                nuevo=caso.union(caso2)
                                lista2.add(nuevo)
                            if len(caso.union(caso2))==iMio:
                                nuevo=caso.union(caso2)
                                finales.add(nuevo)
                                indice[nuevo]=t

                lista2=reduceLista2()
                lista4=set()
                for C in lista2:
                    if C in powers:
                        po2=powers[C]
                    else:
                        po2=set([frozenset(a) for a in mi.powerset(C) if len(a)>0])
                        powers[C]=po2
                    if po2.isdisjoint(todosMCS):
                        lista4.add(C)
                cont2+=1
                lista=lista4
                lista3=[C for C in lista if len(C)<i]

            finales2=finales.copy()
            ########################## He terminado de procesar los soportes ###############
            
            lista3=[C for C in lista if len(C)<i]
            
            misFinales=list(finales)

            lista=set()
            numLista=8*10**6
            r=0
            while r*numLista<len(misFinales):
                maximo=(r+1)*numLista
                if len(misFinales)<maximo:
                    maximo=len(misFinales)
                finales=misFinales[r*numLista:maximo]
                r+=1

                finales=list(finales)

                cpus=cpus1
                k=0
                argus=[]
                num=1000

                while k*num<len(finales):
                    argus.append((k*num,(k+1)*num))
                    k+=1

                with mp.Pool(processes=cpus1) as mp_pool:
                    results=mp_pool.starmap(chequea,argus)

                for result in results:
                    lista.update(result)
            
            del misFinales
            del results
            del finales
            del listaSoportes
            gc.collect()
            
            cont2=1       
            listaL2=list(lista)
            nuevos=set()
            lista2=set()
            
            k=0
            argus=[]
            
            num=miNum1
            cpus=cpus1
            
            while k*num<len(listaL2):
                argus.append((k*num,(k+1)*num+1))
                k+=1

            with mp.Pool(processes=cpus1) as mp_pool:
                results=mp_pool.starmap(analizaSalto,argus)

            for result in results:
                lista2.update(result[0])
                misSoportes.update(result[1])
                nuevos.update(result[1])  
                cont2+=1

            for sop in nuevos:
                if not sop in sopG:
                    sopG[sop]=soporteG(sop,tope)

            listaL2=list(lista2)
            for C in lista2:
                if not C in powers:
                    powers[C]=set([frozenset(a) for a in mi.powerset(C) if len(a)>0])

            lista3=set()
            k=0
            argus=[]
            
            num=1000
            while k*num<len(listaL2):
                argus.append((k*num,(k+1)*num))
                k+=1

            with mp.Pool(processes=cpus1) as mp_pool:
                results=mp_pool.starmap(filtra,argus)

            for result in results:
                lista3.update(result)

            lista2=lista3
            listaL2=list(lista2)
            
            cutsetsT=set()
            
            if i==tope:
                miLista2=lista2

            descartados=set()
            cont=0
            cpus=cpus1
            miLista2=list(lista2)
            listaT=set()

            k=0
            argus=[]
            
            num=200
            while k*num<len(miLista2):
                argus.append((k*num,(k+1)*num))
                k+=1

            with mp.Pool(processes=cpus1) as mp_pool:
                results=mp_pool.starmap(analiza,argus)

            for result in results:
                listaT.update(result[0])
                misSoportes.update(result[1])
                for sop in result[1]:
                    if not sop in sopG:
                        sopG[sop]=soporteG(sop,tope)
                        
            cutsets.append(set(listaT))
            todosMCS.update(listaT)   
        
    sopG=dict()
    powers=dict()
    misSoportes=set()
    return(cutsets)

def analizaSalto(inicio,fin):
    misSoportes2=set()
    descartados=set()
    miLista2=set()
    
    i2=inicio
    while i2 < fin:
        casos=[]
        contCasos=0
        for k in range(i2,fin):
            if k<len(listaL2) and not listaL2[k] in descartados:
                casos.append(listaL2[k])
                contCasos+=1
            if contCasos>=salto:
                break
        i2=k+1

        cont3=0
        if len(casos) > 0:
            with modelo:
                for caso in casos:
                    descartados.add(caso)
                    for g in caso:
                        modelo.genes.get_by_id(g).knock_out()
                #objetivo.bounds=[0,10000]
                sol=modelo.slim_optimize(error_value=-1)
                if sol<=10**-8:
                    for caso in casos:
                        miLista2.add(caso)
                else:
                    sol2=modelo.optimize()
                    sop=frozenset([j for j in range(len(modelo.reactions)) if abs(sol2.fluxes.iloc[j])>10**-8])
                    if not sop in sopG:
                        sopG[sop]=soporteG(sop,tope)
                    misSoportes2.add(sop)
                    for C in set(listaL2).difference(descartados): 
                        sst=time.time()
                        if C in powers:
                            sub=powers[C]
                        else:
                            sub=set([frozenset(a) for a in mi.powerset(C) if len(a)>0])
                            powers[C]=sub
                        if sub.isdisjoint(sopG[sop]):
                            descartados.add(C)
                
    return([miLista2,misSoportes2])



def analiza(inicio,fin):    
    misSoportes2=set()
    descartados=set()
    CSs=set()
    for caso in miLista2[inicio:fin]:
        if not caso in descartados:
            #descartados.add(caso)
            with modelo:
                for g in caso:
                    modelo.genes.get_by_id(g).knock_out()
                #objetivo.bounds=[0,10000]
                sol=modelo.slim_optimize(error_value=-1)
                if sol<=10**-8:
                    CSs.add(caso)
                elif iMio<tope:
                    sol2=modelo.optimize()
                    sop=frozenset([j for j in range(len(modelo.reactions)) if abs(sol2.fluxes.iloc[j])>10**-8])
                    if not sop in sopG:
                        sopG[sop]=soporteG(sop,tope)
                    misSoportes2.add(sop)
                    for C in set(miLista2).difference(descartados): 
                        sst=time.time()
                        sub=set([frozenset(a) for a in mi.powerset(C) if len(a)>0])
                        #if sub.isdisjoint(sopG[sop]):
                        #    descartados.add(C)
                        
    return([CSs,misSoportes2])



def analizaGenes(inicio,fin):
    correctas=set()
    soportes=set()
    for i in range(inicio,fin+1):
        if i<len(listaGenes):
            g=listaGenes[i]
            with modelo:
                modelo.genes.get_by_id(g).knock_out()
                #objetivo.bounds=[0,10000]
                sol=modelo.slim_optimize(error_value=-1)
                if sol<=10**-8:
                    correctas.add(frozenset([g]))
                elif tope>=1:
                    sol2=modelo.optimize(objective_sense="maximize")
                    sop=frozenset([j for j in range(len(modelo.reactions)) if abs(sol2.fluxes.iloc[j])>10**-8])
                    if len(sop)>0:
                        soportes.add(sop)
                    else:
                        correctas.add(frozenset([g]))
    return([correctas,soportes])




def filtra(inicio,fin):    
    CSs=set()
    for caso in listaL2[inicio:fin]:
        vale=True
        po=powers[caso]
        for sop in nuevos:
            if po.isdisjoint(sopG[sop]):
                vale=False
                break
        if vale:
            CSs.add(caso)
                        
    return(CSs)


def createGMCS(_model,_tope,_solver="gurobi",_cpus0=64,_cpus1=64,_salto=12,_objective=""):
    global modelo, tope, listaGenes, misSoportes, salto, asociados, sopG, powers
    global cpus0, cpus1, objetivo, miNum0, miNum1
    
    cobra.Configuration().solver=_solver 
    miNum0=30
    miNum1=2000
    modelo=_model
    tope=_tope
    listaGenes=list(modelo.genes)
    
    misSoportes=set()
    
    salto=_salto
    cpus0=_cpus0
    cpus1=_cpus1
        
    if len(modelo.reactions)<=1000:
        salto=4

    asociados,reacciones,eliminar=asociar(modelo,tope,listaGenes)
    
    if _objective=="":
        if len(cobra.util.solver.linear_reaction_coefficients(modelo))==1:
             for r in cobra.util.solver.linear_reaction_coefficients(modelo):
                objetivo=r
    else:
        objetivo=modelo.reactions.get_by_id(_objective)
    
    
    sopG=dict()
    powers=dict()
    sol=modelo.slim_optimize(error_value=-1)
    if sol==-1:
        print("Fallo")
        CS1=[set(),set()]
    else:
        CS1=CSV2(modelo,tope)
    return(CS1)