#Lanza programa FreeFem++ ElasticityDamage3D.edp con diferentes datos


cat > Defaults.idp <<END
VersionEta      =0; 
UsoEpsilonUpunto=1; 
VersionDivEta   =1; 
END

FreeFem++ ElasticityDamage3D.edp > 011_salida.log &
echo pausa
ping 127.0.0.1 -n 2 > nul

cat > Defaults.idp <<END
VersionEta      =1; 
UsoEpsilonUpunto=1; 
VersionDivEta   =1; 
END

FreeFem++ ElasticityDamage3D.edp > 111_salida.log &
echo pausa
ping 127.0.0.1 -n 2 > nul



cat > Defaults.idp <<END
VersionEta      =0; 
UsoEpsilonUpunto=0; 
VersionDivEta   =1; 
END

FreeFem++ ElasticityDamage3D.edp > 001_salida.log &
echo pausa
ping 127.0.0.1 -n 2 > nul

cat > Defaults.idp <<END
VersionEta      =1; 
UsoEpsilonUpunto=0; 
VersionDivEta   =1; 
END

FreeFem++ ElasticityDamage3D.edp > 101_salida.log &
echo pausa
ping 127.0.0.1 -n 2 > nul
