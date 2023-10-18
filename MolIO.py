#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import os
import sys
import tempfile
from pathlib import Path

import cmder
from loguru import logger

LOG_LEVEL = 'ERROR'
logger.remove()
logger.add(sys.stdout, colorize=True, format='<level>[{time:YYYY-MM-DD HH:mm:ss}] {message}</level>', level=LOG_LEVEL)


class SDF:
    def __init__(self, s):
        self.s = s
        self.title, self.mol, self.properties, self.properties_dict = self.parse(s)

    @staticmethod
    def parse(s):
        if s:
            lines = s.splitlines(keepends=True)
            title = lines[0].rstrip()

            mol, i = [], 0
            for i, line in enumerate(lines[1:]):
                mol.append(line)
                if line.strip() == 'M  END':
                    break
            mol = ''.join(mol)

            properties = lines[i + 2:]

            properties_dict, idx = {}, []
            for i, p in enumerate(properties):
                if p.startswith('>'):
                    properties_dict[p.split(maxsplit=1)[1].rstrip()] = ''
                    idx.append(i)
            idx.append(-1)

            for i, k in enumerate(properties_dict.keys()):
                properties_dict[k] = ''.join(properties[idx[i] + 1:idx[i + 1]])

            properties = ''.join(properties[:-1])
        else:
            title, mol, properties, properties_dict = '', '', '', {}
        return title, mol, properties, properties_dict

    def sdf(self, output=None, properties=None, title=''):
        if properties:
            properties = self.properties if isinstance(properties, bool) else properties
            if isinstance(properties, dict):
                properties = '\n'.join([f'> {k}\n{v}' for k, v in properties.items()])
            else:
                logger.error('Invalid properties, only accept a boolean value or a dictionary')
                sys.exit(1)
            s = f'{title or self.title}\n{self.mol}\n{properties}\n$$$$\n'
        else:
            s = f'{title or self.title}\n{self.mol}\n$$$$\n'

        if output:
            with open(output, 'w') as o:
                o.write(s)
            return output
        else:
            return s

    def __str__(self):
        return self.s
    
    
class DLG:
    def __init__(self, s):
        self.s = s
        self.user, self.remark, self.atom = self.parse(s)

    @staticmethod
    def parse(s):
        user, remark, atom = [], [], []
        prefix, lines = 'DOCKED: ', s.splitlines(keepends=True)
        for i, line in enumerate(lines):
            if line.startswith(prefix):
                line = line.removeprefix(prefix)
                if line.startswith('USER'):
                    user.append(line)
                elif line.startswith('REMARK'):
                    remark.append(line)
                elif line.startswith('ATOM'):
                    atom.append(line)
        return [''.join(x) for x in (user, remark, atom)]
    
    def pdbqt(self, output=None, properties=None, title=''):
        if properties:
            properties = self.remark.splitlines() if isinstance(properties, bool) else properties
            if isinstance(properties, list):
                properties = [p if p.startswith('REMARK') else f'REMARK  {p}' for p in properties]
            else:
                logger.warning(f'Invalid properties {properties}, properties not added to output')
            if title:
                for p in properties:
                    if 'REMARK  Name == ' in p:
                        break
                else:
                    properties.insert(0, f'REMARK  Name = {title}')
            s = '\n'.join([self.user, '\n'.join(properties), self.atom])
        else:
            s = '\n'.join(line.removeprefix('DOCKED: ') for line in self.s.splitlines()
                          if line.startswith('DOCKED: ATOM'))
            if title:
                if 'REMARK  Name = ' not in s:
                    s = f'REMARK  Name = {title}\n{s}'
        if output:
            with open(output, 'w') as o:
                o.write(s)
                return output
        else:
            return s

    def sdf(self, output=None, properties=None, title=''):
        pdbqt, _sdf = tempfile.mktemp(suffix='.pdbqt'), output or tempfile.mktemp(suffix='.sdf')
        self.pdbqt(output=pdbqt, properties=properties, title=title)
        p = cmder.run(f'obabel {pdbqt} -o sdf -O {_sdf}', exit_on_error=False)
        if p.retruncode:
            os.unlink(pdbqt)
            os.unlink(_sdf)
            logger.error('Failed to write molecule to SDF file')
            return ''

        os.unlink(pdbqt)
        if output:
            return output
        else:
            with open(output) as f:
                s = f.read()
            os.unlink(output)
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


def parse(path):
    if path.endswith('.sdf.gz') or path.endswith('.sdf'):
        return parse_sdf(path)
    elif path.endswith('.dlg.gz') or path.endswith('.dlg'):
        return parse_dlg(path)
    else:
        raise TypeError(f'Unsupported file {path}')


def write(ss, output, properties=False):
    with open(output, 'w') as o:
        o.writelines(f'{s.sdf(properties=properties) if isinstance(s, SDF) else s}\n$$$$\n' for s in ss)
    return output


if __name__ == '__main__':
    for item in parse('tests/autodock.dlg'):
        s = item.pdbqt(title='000L-5787', properties=['score = 5.03e+07'], output='000L-5787.pdbqt')
        print(s)
        break
