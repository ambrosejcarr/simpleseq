#!/usr/local/bin/python3

import regex as re


_pattern = re.compile(b'(.{8,11}?)(GAGTGATTGCTTGTGACGCCTT){s<=1}(.{8})(.{6})(.*)')


def in_drop(g, b):
    try:
        cell1, spacer, cell2, rmt, poly_t = re.match(
            _pattern, b.sequence).groups()
        cell = cell1 + cell2
    except AttributeError:
        cell, rmt, poly_t = b'', b'', b''
    return g.add_annotation((b'', cell, rmt, poly_t))


def drop_seq(g, b):
    cell = b.sequence[:12]
    rmt = b.sequence[12:20]
    return g.add_annotation((b'', cell, rmt, b''))


def mars1_seq(g, b=None):
    *name_fields, pool, cell, rmt = g.name[:-1].split(b':')
    g.name = b':'.join(name_fields) + b';' + b':'.join((pool, cell, rmt, b'')) + b'\n'
    return g


def mars2_seq(g, b=None):
    *name_fields, pool, cell, rmt = g.name[:-1].split(b'-')
    g.name = b'-'.join(name_fields) + b';' + b':'.join((pool, cell, rmt, b'')) + b'\n'
    return g
