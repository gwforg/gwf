.PHONY: init test lint coverage docs

init:
	pip install pipenv --upgrade
	pipenv install --dev --skip-lock

test:
	pipenv run coverage run --source gwf setup.py test

lint:
	pipenv run flake8 src/gwf

coverage:
	pipenv run coverage report

docs:
	pipenv run $(MAKE) -C docs html

clean:
	rm -rf docs/_build .gwfconf.json
