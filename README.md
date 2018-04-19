# Option Calculator Website

This repository contains sources of Django application that powers Option Calculator.

## What's in it?

It's a simple calculator that contains 7 models:

- European call/put option
- Implied volatility calculator
- American call/put option
- Geometric Asian option
- Arithmetic Asian option
- Geometric basket option
- Arithmetic basket option

# Run the website

Please note that we use Python 2.7, so make sure that you use correct version when running commands below.

## Setting up a development environment

First, clone the repository:
Step into newly created `djangogirls` directory:

    cd ass3_website

Create a new virtual environment if needed. Then, install all the required dependencies:

    pip install -r requirements.txt

Run the migration

    python manage.py makemigrations
    python manage.py migrate

Run your local server:

    python manage.py runserver


