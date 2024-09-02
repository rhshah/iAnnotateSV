# options
.ONESHELL:


.PHONY: deps-install
deps-install:  ## install dependencies
	pip install poetry
	poetry install 

requirements.txt: poetry.lock
	poetry export --format requirements.txt --output requirements.txt --without-hashes

requirements-dev.txt: poetry.lock
	poetry export --with dev --format requirements.txt --output requirements-dev.txt --without-hashes