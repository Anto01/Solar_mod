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

I = 
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