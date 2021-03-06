.DEFAULT_GOAL := all

force-build:
	@docker pull ubuntu:latest
	@docker build -t hivdb/gaggp41-runtime:latest --force-rm --no-cache .

build:
	@docker build -t hivdb/gaggp41-runtime:latest .

shell:
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime /bin/bash

msa:
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/process-msa.sh

basic: msa
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-basic.sh

nj:
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-nj.sh

naive:
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime python3 scripts/analyze-naive-studies.py

naive-dist:
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-naive-dist.sh

pairwise: basic
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-pairwise.sh

# relax: basic nj
# 	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-relax.sh

fel: basic nj
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-fel.sh

meds: basic nj
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-meds.sh

final:
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/run-final.sh

report:
	@docker run -it --rm --volume `pwd`:/app hivdb/gaggp41-runtime scripts/generate-report.sh

all: basic pairwise nj fel meds final report
.PHONY: report
