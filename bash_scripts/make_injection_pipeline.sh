$SOURCE_INJECTION_DIR/bin/make_injection_pipeline \
	-t deepCoadd \
	-r $DRP_PIPE_DIR/pipelines/HSC/DRP-RC2.yaml \
	-i $SOURCE_INJECTION_DIR/pipelines/inject_coadd.yaml \
	-f DRP-RC2+injection_by_peter.yaml
