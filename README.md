ElasticityDamage3D.edp
======================
Programa principal en FreeFem++

Defaults.idp 
============
Se usa para cambiar algunos datos para la ejecucion.
Ejemplo:
VersionEta      =1;
UsoEpsilonUpunto=0;
VersionDivEta   =1;

LanzaCalculos.sh
================
Script para lanzar varios cálculos a la vez (cada cálculo usa solo 1 core)
La idea de este Script, es realizar 2 etapas:
	modificar el archivo Defaults.idp,
	Lanzar FreeFem++ con esta modificación.
	
Cada vez que FreeFem se lanza, genera un log y archivos de resultados 
que tienen un prefijo para identificarlo. Para más informacion, ver
la linea 20 de ElasticityDamage3D.edp, donde se define la variable prefix
como:
	string prefix=  ""+VersionEta+  UsoEpsilonUpunto+VersionDivEta+"_";