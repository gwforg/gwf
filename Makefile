init:
	pip install -r requirements.txt
	pip install .
test:
	coverage run --source gwf -m unittest discover tests/
lint:
	flake8
