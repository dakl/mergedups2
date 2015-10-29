# mergedups2

## Install

git clone https://github.com/dakl/mergedups2
pip install ./mergedups2

## Test

~~~bash
py.test --cov=mergedups --cov=tests/test_pileup.py --cov tests/test_readmerger.py --cov-report=html -v
~~~

## Run

~~~bash
mergedups2 -I test.bam -x metrics.txt -O /dev/stdout | samtools sort - testout-sort
~~~

## Debug

~~~bash
mergedups2 --loglevel DEBUG -I test.bam --reads_between_logs 1 -x testout.metrics -O /dev/stdout | samtools sort - testout-sort
~~~