# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 20:17:45 2016

@author: Antoine
"""

from sys import *
from solar_mod import *
from numpy import *
from collections import *


# Question #1

temps = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
I = [0,0,0,0,0,27,99,280,432,543,636,802,796,791,609,396,339,226,83,16,0,0,0]
Ibn = [0,0,0,0,0,98,110,411,499,423,387,579,521,519,578,348,142,175,285,86,0,0,0,0]
thez = 60
rhog = 0.4
lon = -73.0 -35.0/60.0 # longitude  Montreal (73 deg 35 min ouest)
Lst = -75 # longitude  méridien HNE -5
phi  = 45.0 +30.0/60.0  # latitude  Montreal (45 deg 30 min nord)  

      
st  = namedtuple('temps',['hr','min','jour'])


# 
n =jour_mois_jour_annee(12,'juin')
print ('n=',n)

del_h = 1.0  # heure avancée

### Mettre ça en boucle
st.hr = 12.0
st.min = 0.0
st.jour = n
solt = heure_solaire(lon,Lst,del_h,st)



print ('heure solaire = '),solt.hr
print ('minute solaire = '),solt.min
print ('jour  solaire = '),solt.jour


Ib = irradiation_horaire_tilted(Ibn,thez)



ome1 = -3.0*15.0
ome2 = ome1+15.0

sol_t = heure_solaire(lon,Lst,del_h,st=1)

angle_horaire(sol_t)



Io = irradiation_extraterrestre_horaire(n,phi,ome1,ome2)
kt = indice_clarete_horaire(I,Io)

# Calcul pour le 20 février
#
jour = 20
mois = 'fev'
n = jour_mois_jour_annee(jour,mois);
ome1 = -3.0*15.0
ome2 = ome1+15.0
Io = irradiation_extraterrestre_horaire(n,phi,ome1,ome2)
kt = I/Io