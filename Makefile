.DEFAULT_GOAL := all

build:
	@docker build -t gaggp41-runtime .

shell:
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime /bin/bash

conda:
	@docker run -it --rm --volume `pwd`:/app continuumio/miniconda /bin/bash

basic:
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime scripts/run-basic.sh

nj:
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime scripts/run-nj.sh

pairwise: basic
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime scripts/run-pairwise.sh

relax: basic nj
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime scripts/run-relax.sh

fel: basic nj
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime scripts/run-fel.sh

meds: basic nj
	@docker run -it --rm --volume `pwd`:/app gaggp41-runtime scripts/run-meds.sh

all: basic pairwise nj relax fel meds
