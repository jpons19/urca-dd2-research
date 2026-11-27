.PHONY: test regenerate-t9 regenerate-t8 multi-eos all

test:
	python calculations/test_final.py
	python calculations/test_T8_reference.py
	python -m unittest discover -s tests -v

regenerate-t9:
	python calculations/build_final.py

multi-eos:
	python -m calculations.multi_eos.run_compare --points 120

regenerate-t8:
	python calculations/run_dd2_T8_reference.py

all: regenerate-t9 regenerate-t8 test
