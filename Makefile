.DEFAULT_GOAL := all

build:
	@docker build -t gaggp41-runtime .

shell:
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime /bin/bash

all:
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime scripts/run.sh
