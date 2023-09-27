DATA_COLL=HSC/runs/RC2/w_2023_32/DM-40356
REPO=/repo/main

names='27 29 30 31 33'
for name in $names
do
    INJECT_COLL=u/pferguso/full_dwarf_test/2_lim_$name
    OUTPUT=u/pferguso/full_dwarf_test/2_lim_$name/step3
    echo $INJECT_COLL
    pipetask --long-log \
    run --register-dataset-types \
    --instrument lsst.obs.subaru.HyperSuprimeCam \
    -b $REPO \
    -i $DATA_COLL,$INJECT_COLL \
    -o $OUTPUT \
    -p pipelines/DRP-RC2+injection_by_peter.yaml#step3 \
    -d "instrument='HSC' AND skymap='hsc_rings_v1' AND
    tract=9615 AND patch=3 AND (band = 'g' OR band = 'r' OR band = 'i') " \
    -j 6
done

