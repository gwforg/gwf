.PHONY: init test lint coverage docs

init:
	pip install -r requirements.txt
	pip install -e .
test:
	coverage run --source gwf setup.py test
lint:
	flake8
coverage:
	coverage report -m
docs:
	$(MAKE) -C docs html
