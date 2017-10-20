#Lanza programa FreeFem++ ElasticityDamage3D.edp con diferentes datos

version="new"
casos=(
#   "001"  # (0)eta=0       , (0)E(u),    (1)    div(u)
#   "002"  # (0)eta=0       , (0)E(u),    (2)abs(div(u))
   "003"  # (0)eta=0       , (0)E(u),    (3)max(div(u),0)
#   "011"  # (0)eta=0       , (1)E(udot), (1)    div(u)
#   "012"  # (0)eta=0       , (1)E(udot), (2)abs(div(u))
   "013"  # (0)eta=0       , (1)E(udot), (3)max(div(u),0)
#   "101"  # (1)eta=articulo, (0)E(u),    (1)    div(u)
#   "102"  # (1)eta=articulo, (0)E(u),    (2)abs(div(u))
   "103"  # (1)eta=articulo, (0)E(u),    (3)max(div(u),0)
#   "111"  # (1)eta=articulo, (1)E(udot), (1)    div(u)
#   "112"  # (1)eta=articulo, (1)E(udot), (2)abs(div(u))
   "113"  # (1)eta=articulo, (1)E(udot), (3)max(div(u),0)
   "201"  # (2)eta=MGrande , (0)E(u),    (1)    div(u)
#   "202"  # (2)eta=MGrande , (0)E(u),    (2)abs(div(u))
   "203"  # (2)eta=MGrande , (0)E(u),    (3)max(div(u),0)
   "211"  # (2)eta=MGrande , (1)E(udot), (1)    div(u)
#   "212"  # (2)eta=MGrande , (1)E(udot), (2)abs(div(u))
   "213"  # (2)eta=MGrande , (1)E(udot), (3)max(div(u),0)
)

case $version in
    "new")
        for i in "${casos[@]}"
        do
            echo $i
            echo $i|
                sed -e 's/^\(.\)\(.\)\(.\)/VersionEta =\1;UsoEpsilonUpunto =\2;VersionDivEta =\3;/g'>Defaults.idp
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