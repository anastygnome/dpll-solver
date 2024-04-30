# This file is part of the anastygnome repository
# Copyright (c) 2022-23 Tom Domenge.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_file(filename):
    def read_clause(file, nb_lit):
        buffer = ""
        for line in file:
            if not line.strip():
                continue
            buffer += line.strip()
            if line.strip().endswith("0"):
                current_clause = tuple(map(int, buffer.rstrip("0").split()))
                if any(abs(literal) > nb_lit for literal in current_clause):
                    raise ValueError("Invalid clause")
                yield current_clause
                buffer = ""  # Reset the buffer after processing a clause
            else:
                buffer += " "

    cnf = set()
    nb_clauses = 0
    nb_lit = 0
    with filename.open() as f:
        nb_clauses, nb_lit = map(int, f.readline().strip().split())
        print(
            f"Reading {nb_clauses} clause(s) with at most {nb_lit} literal(s) each from {filename}"
        )
        for clause in read_clause(f, nb_lit):
            cnf.add(clause)
            if len(cnf) >= nb_clauses:
                break

    if len(cnf) < nb_clauses:
        print(
            f"Warning, {len(cnf)} clauses in file is less than ({nb_clauses}) expected, the CNF may be truncated."
        )
    return cnf, nb_clauses, nb_lit


def pure_literal(cnf, symbols, model):
    if not cnf:
        return None, None
    polarity = {}
    new_symb = set() if not symbols else symbols.copy()
    pure_literals = set() if not symbols else symbols.copy()

    for clause in cnf:
        for lit in clause:
            abs_lit = abs(lit)
            if abs_lit not in polarity:
                polarity[abs_lit] = lit > 0
                if not symbols:
                    new_symb.add(abs_lit)
                    pure_literals.add(abs_lit)
            elif polarity[abs_lit] != (lit > 0):
                pure_literals.discard(abs_lit)

    new_symb -= pure_literals

    for lit in pure_literals:
        model[lit] = polarity[lit]
    return {
        clause for clause in cnf if not any(abs(lit) in pure_literals for lit in clause)
    }, new_symb


def remove_unit_clauses(clauses, model):
    """Remove unit clauses and update the model."""
    for clause in clauses.copy():
        unit_literal = None
        for literal in clause:
            if abs(literal) not in model:
                if (
                    unit_literal is not None
                ):  # More than one unassigned literal in the clause
                    break
                unit_literal = literal
        else:  # No break occurred, meaning only one unassigned literal found in the clause
            if unit_literal is not None:
                model[abs(unit_literal)] = unit_literal > 0
                clauses.remove(clause)
                return True
    return False


def DPLL(clauses, symbols=None, model={}):
    if not clauses:
        return (True, None)

    # Check for all clauses satisfied by the model (this actually computes the truth value of the formula)
    if model and all(
        any(
            (
                litteral in model
                and model[litteral]
                or -litteral in model
                and not model[-litteral]
                for litteral in c
            )
        )
        for c in clauses
    ):
        return (True, model)
    elif symbols is not None and len(model) == len(symbols):
        return (False, model)
    if remove_unit_clauses(clauses, model):
        return DPLL(clauses, symbols, model)

    clauses, symbols = pure_literal(clauses, symbols, model)
    if not clauses:
        return (True, model)

    # Choose an unassigned symbol and branch
    for s in symbols - model.keys():
        model[s] = True
        result = DPLL(clauses, symbols, model)
        if result:
            return (True, model)

        model[s] = False
        return DPLL(clauses, symbols, model)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check if the CNF described in a file is satisfiable."
    )
    parser.add_argument("filename", type=Path, help="The file to read the CNF from.")
    args = parser.parse_args()

    cnf, nb_clauses, nb_lit = parse_file(args.filename)
    result = DPLL(cnf)

    if result:
        print(f"Satisfiable with model: {result[1]}")
    else:
        print("Unsatisfiable")
