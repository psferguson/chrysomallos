DATA_COLL=HSC/runs/RC2/w_2023_32/DM-40356
REPO=/repo/main

bands='g'
names='27 34 28 29 30 31 32 33'
for band in bands
do
    for name in $names
        do
            INJECT_COLL=u/pferguso/maglim_16_test/maglim_${name}_round_1
            OUTPUT=u/pferguso/maglim_16_test/maglim_${name}_round_1/injection_step
            echo $INJECT_COLL
            pipetask --long-log \
            run --register-dataset-types \
            --instrument lsst.obs.subaru.HyperSuprimeCam \
            -b $REPO \
            -i $DATA_COLL,$INJECT_COLL \
            -o $OUTPUT \
            -p pipelines/DRP-RC2+injection_by_peter.yaml#inject_coadd \
            -d "instrument='HSC' AND skymap='hsc_rings_v1' AND
            tract=9615 AND patch=3 AND (band = '$band') " 
            
        done
done
