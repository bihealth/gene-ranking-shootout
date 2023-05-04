.PHONY: default
default:

.PHONY: install-dev
install-dev:
	pip install -U pip setuptools
	pip install -U "black==22.3.0" "isort>=5.0,<6.0" "flake8>=5.0,<6.0" mypy pytest pytest-cov \
		types-requests types-tqdm
	pip install -e .

.PHONY: black
black:
	black -l 100 .

.PHONY: black-check
black-check:
	black -l 100 --check .

.PHONY: isort
isort:
	isort --force-sort-within-sections --profile=black .

.PHONY: isort-check
isort-check:
	isort --force-sort-within-sections --profile=black --check .

.PHONY: flake8
flake8:
	flake8

.PHONY: mypy
mypy: export MYPYPATH=stubs
mypy:
	mypy gene_ranking_shootout

.PHONY: lint
lint: flake8 isort-check black-check mypy

.PHONY: pytest
pytest:
	pytest .
