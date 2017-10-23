#Lanza programa FreeFem++ ElasticityDamage3D.edp con diferentes datos

version="new"
casos=(
   "0020"  # (0)eta=0       , (0)E(u),    (1)    div(u)    ,(0) alpha=0
#   "0011"  # (0)eta=0       , (0)E(u),    (1)    div(u)    ,(1) alpha evoluciona
#   "0011"  # (0)eta=0       , (0)E(u),    (1)    div(u)
#   "0021"  # (0)eta=0       , (0)E(u),    (2)abs(div(u))
#   "0031"  # (0)eta=0       , (0)E(u),    (3)max(div(u),0)
#   "0111"  # (0)eta=0       , (1)E(udot), (1)    div(u)
#   "0121"  # (0)eta=0       , (1)E(udot), (2)abs(div(u))
#   "0131"  # (0)eta=0       , (1)E(udot), (3)max(div(u),0)
#   "1011"  # (1)eta=articulo, (0)E(u),    (1)    div(u)
#   "1021"  # (1)eta=articulo, (0)E(u),    (2)abs(div(u))
#   "1031"  # (1)eta=articulo, (0)E(u),    (3)max(div(u),0)
#   "1111"  # (1)eta=articulo, (1)E(udot), (1)    div(u)
#   "1121"  # (1)eta=articulo, (1)E(udot), (2)abs(div(u))
#   "1131"  # (1)eta=articulo, (1)E(udot), (3)max(div(u),0)
#   "2011"  # (2)eta=MGrande , (0)E(u),    (1)    div(u)
#   "2021"  # (2)eta=MGrande , (0)E(u),    (2)abs(div(u))
#   "2031"  # (2)eta=MGrande , (0)E(u),    (3)max(div(u),0)
#   "2111"  # (2)eta=MGrande , (1)E(udot), (1)    div(u)
#   "2121"  # (2)eta=MGrande , (1)E(udot), (2)abs(div(u))
#   "2131"  # (2)eta=MGrande , (1)E(udot), (3)max(div(u),0)
)

case $version in
    "new")
        for i in "${casos[@]}"
        do
            echo Calculando caso $i
            echo $i|
                sed -e \
                's/^\(.\)\(.\)\(.\)\(.\)/VersionEta =\1;UsoEpsilonUpunto =\2;VersionDivEta =\3;UsoDalphaDt=\4;/g'>Defaults.idp
            cat Defaults.idp
            echo 'FreeFem++ ElasticityDamage3D.edp > '$i'_salida.log &'
            FreeFem++ ElasticityDamage3D.edp > ${i}_salida.log &
            echo pausa
            ping 127.0.0.1 -n 2 > nul
        done
        ;;

    "old")  # Era lo mismo que "new", pero repitiendo las lineas una a una....
        cat > Defaults.idp <<END
        VersionEta      =0; UsoEpsilonUpunto=1; VersionDivEta   =2; 
END

        FreeFem++ ElasticityDamage3D.edp > 012_salida.log &
        echo pausa
        ping 127.0.0.1 -n 2 > nul

        cat > Defaults.idp <<END
        VersionEta      =1; UsoEpsilonUpunto=1; VersionDivEta   =2; 
END

        FreeFem++ ElasticityDamage3D.edp > 112_salida.log &
        echo pausa
        ping 127.0.0.1 -n 2 > nul

        ;;


esac