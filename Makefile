.PHONY: init test lint coverage docs

init:
	pip install -r requirements.txt
	pip install -e .

test:
	coverage run --source src/gwf setup.py test
	coverage xml

lint:
	flake8

coverage:
	coverage report

docs:
	$(MAKE) -C docs html
