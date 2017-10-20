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