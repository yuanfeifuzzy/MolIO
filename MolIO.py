#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A simple module for convenient molecule I/O
"""

import os
import sys

import cmder
from loguru import logger

logger.add(sys.stdout, colorize=True, format='<level>[{time:YYYY-MM-DD HH:mm:ss}] {message}</level>', level='INFO')


class SDF:
    def __init__(self, s):
        self.s = s
        self.score, self.title, self.mol = self.parse(s)

    @staticmethod
    def parse(s):
        if s:
            lines = s.splitlines(keepends=True)
            title = lines[0].rstrip()

            mol, i, score = [], 0, None
            for i, line in enumerate(lines[1:]):
                mol.append(line)
                if line.strip() == 'M  END':
                    break
            mol = ''.join(mol)

            for line in lines[i + 2:]:
                if line.startswith('ENERGY') and 'LOWER_BOUND' in line:
                    try:
                        score = line.split('=')[1].lstrip().split()[0]
                    except Exception as e:
                        logger.error(f'Failed to get docking score from {line.strip()} due to {e}')
                        score = None
        else:
            score, title, mol = None, '', ''
        return score, title, mol

    def sdf(self, output=None, title=''):
        score = '' if self.score is None else f'> <score>\n{self.score}'
        if self.mol:
            s = f'{title or self.title}\n{self.mol}'
            if score:
                s = f'{s}\n{score}'
            s = f'{s}\n$$$$\n'
        else:
            s = ''

        if output:
            if s:
                with open(output, 'w') as o:
                    o.write(s)
                return output
            else:
                return ''
        else:
            return s

    def __str__(self):
        return self.s


class DLG:
    def __init__(self, s):
        self.s = s
        self.score, self.title, self.atom = self.parse(s)

    @staticmethod
    def parse(s):
        score, atom, title = None, [], ''
        prefix, lines = 'DOCKED: ', s.splitlines(keepends=True)
        for i, line in enumerate(lines):
            if line.startswith(prefix):
                line = line.removeprefix(prefix)
                if line.strip().split()[0] in ('ATOM', 'ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH'):
                    atom.append(line)
                elif 'Estimated Free Energy of Binding' in line:
                    score = line.split('=')[1].lstrip().split()[0]
                elif 'REMARK Name =' in line:
                    title = line.strip().split(' = ')[1]

        return score, title, ''.join(atom)

    def pdbqt(self, output=None, title=''):
        title = title or self.title
        s = [f'REMARK  Name = {title}' if title else '', '' if self.score is None else f'score {self.score}', self.atom]
        s = '\n'.join([x for x in s if x])
        if output:
            with open(output, 'w') as o:
                o.write(s)
                return output
        else:
            return s

    def sdf(self, output=None, title=''):
        pdbqt, _sdf = tempfile.mktemp(suffix='.pdbqt'), output or tempfile.mktemp(suffix='.sdf')
        self.pdbqt(output=pdbqt, title=title)
        p = cmder.run(f'obabel {pdbqt} -b -p -h -n -o sdf -O {_sdf}', exit_on_error=False)
        if p.returncode:
            os.unlink(pdbqt)
            os.unlink(_sdf)
            logger.error('Failed to write molecule to SDF file')
            return ''

        os.unlink(pdbqt)
        if output:
            return _sdf
        else:
            with open(_sdf) as f:
                s = f.read()
            os.unlink(_sdf)
            return s

    def __str__(self):
        return self.s


def parse_sdf(sdf):
    path = str(sdf)
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        lines = []
        for i, line in enumerate(f):
            lines.append(line)
            if line.strip() == '$$$$':
                yield SDF(''.join(lines))
                lines = []
                continue
        yield SDF(''.join(lines))


def parse_dlg(dlg):
    path = str(dlg)
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        lines = []
        for line in f:
            if line.startswith('DOCKED: MODEL'):
                break
        for line in f:
            if line.startswith('DOCKED: ENDMDL'):
                yield DLG(''.join(lines))
                lines = []
                continue
            lines.append(line)

        yield DLG(''.join(lines))


def write(records, output=''):
    records = (record.sdf() for record in records)
    s = ''.join(record for record in records if record)
    if output:
        with open(output, 'w') as o:
            o.write(s)
        return output
    else:
        return s


def split_sdf(sdf, prefix, suffix='.sdf', files=0, records=0):
    names = []
    if files:
        n = sum(1 for _ in parse_sdf(sdf))
        num, remains = divmod(n, files)
        if remains:
            num += 1

        idx, end, items = 1, num, []
        for i, item in enumerate(parse_sdf(sdf), 1):
            items.append(item)
            if i == end:
                end += num
                name = write(items, output=f'{prefix}.{idx}.{suffix}')
                names.append(name)
                idx += 1
                items = []
        if items[:-1]:
            name = write(items[:-1], output=f'{prefix}.{idx}.{suffix}')
            names.append(name)
    elif records:
        idx, n, items = 0, records - 1, []
        for i, item in enumerate(parse_sdf(sdf), 0):
            items.append(item)
            idx, remains = divmod(i, records)

            if remains == n:
                name = write(items, output=f'{prefix}.{idx + 1}.{suffix}')
                names.append(name)
                items = []
        if items[:-1]:
            name = write(items[:-1], output=f'{prefix}.{idx + 1}.{suffix}')
            names.append(name)
    else:
        idx = 1
        name = write(parse_sdf(sdf), output=f'{prefix}.{idx}.{suffix}')
        names.append(name)
    return names


def clean_sdf(sdf, output=None):
    items = [sdf] if isinstance(sdf, SDF) else parse_sdf(sdf)
    s = ''.join(s.sdf() for s in items)
    if output is None:
        return s
    else:
        if output:
            with open(output, 'w') as o:
                o.write(s)
            return output
        else:
            if isinstance(sdf, str):
                output = str(Path(sdf).with_suffix('.clean.sdf'))
                with open(output, 'w') as o:
                    o.write(s)
                return output
            else:
                return s


def dlg2sdf(dlg, sdf=None):
    s = ''.join(f'{item.sdf()}' for item in parse_dlg(dlg))
    pdbqt = tempfile.mktemp(suffix='.pdbqt')
    with open(pdbqt, 'w') as o:
        o.write(s)

    if sdf is None:
        _sdf = tempfile.mktemp(suffix='.sdf')
    else:
        _sdf = sdf or str(Path(dlg).with_suffix('.sdf'))

    p = cmder.run(f'obabel {pdbqt} -b -p -h -n -o sdf -O {_sdf}', exit_on_error=False)
    if p.returncode:
        os.unlink(pdbqt)
        os.unlink(_sdf)
        logger.error('Failed to write molecule to SDF file')
        return ''
    else:
        if sdf is None:
            with open(_sdf) as f:
                s = f.read()
            os.unlink(_sdf)
            return s
        else:
            return _sdf
