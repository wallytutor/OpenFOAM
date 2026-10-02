# -*- coding: utf-8 -*-

import numpy as np

from argparse import ArgumentParser

import methane_air as ma


def main():
    """ CLI workflow for case setup generation. """
    parser = ArgumentParser(
        description = "Prepare case conditions."
    )
    parser.add_argument(
        "--temperature",
        type    = float,
        default = 25.0,
        help    = "Inlet temperature in degrees Celsius."
    )
    args = parser.parse_args()
    ma.load_setup_shared(air_temp=args.temperature, save=True)


if __name__ == "__main__":
    main()
