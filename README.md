ElasticityDamage3D.edp
======================
Programa principal en FreeFem++

Defaults.idp 
============
Se usa para cambiar algunos datos para la ejecucion.<br>
Ejemplo:<br>
<pre>
VersionEta      =1;
UsoEpsilonUpunto=0;
VersionDivEta   =1;
</pre>

LanzaCalculos.sh
================
Script para lanzar varios cálculos a la vez (cada cálculo usa solo 1 core).
La idea de este Script, es realizar 2 etapas:<br>
<li>	modificar el archivo Defaults.idp, 
</li><li>	Lanzar FreeFem++ con esta modificación.
</li>
	
Cada vez que FreeFem se lanza, genera un log y archivos de resultados 
que tienen un prefijo para identificarlo. Para más informacion, ver
la linea 20 de <tt>ElasticityDamage3D.edp</tt>, donde se define la variable prefix
como:<pre>
	string prefix=  ""+VersionEta+  UsoEpsilonUpunto+VersionDivEta+"_";</pre>
    
    Este programa está mas operativo (pero menos comprensible!!!)
    
    Basta con definir un arreglo de prefijos a calcular
y se lanza un for que hace todos los cálculos:

<pre>
version="new"
casos=(
#   "001"  # (0)eta=0       , (0)E(u),    (1)    div(u)
#   "002"  # (0)eta=0       , (0)E(u),    (2)abs(div(u))
    "003"  # (0)eta=0       , (0)E(u),    (3)max(div(u),0)
#   "011"  # (0)eta=0       , (1)E(udot), (1)    div(u)
#   "012"  # (0)eta=0       , (1)E(udot), (2)abs(div(u))
... etc...
</pre>

El for solo ejecuta las lineas no comentadas del arreglo "casos"