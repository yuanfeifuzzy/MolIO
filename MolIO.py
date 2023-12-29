#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A simple module for convenient molecule I/O
"""

import os
import gzip
import random
import tempfile
from pathlib import Path

import cmder
import vstool

logger = vstool.setup_logger()


def _opener(path):
    return gzip.open if str(path).endswith('.gz') else open


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

            for j, line in enumerate(lines[i + 2:]):
                if line.startswith('ENERGY') and 'LOWER_BOUND' in line:
                    try:
                        score = float(line.split('=')[1].lstrip().split()[0])
                    except Exception as e:
                        logger.error(f'Failed to get docking score from {line.strip()} due to {e}')
                        score = None
                elif line.startswith('>') and '<score>' in line:
                    line = lines[i + 2:][j + 1]
                    try:
                        score = float(line.strip())
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
                s = f'{s}\n{score}\n'
            s = f'{s}\n$$$$\n'
        else:
            s = '\n'

        if output:
            if s:
                output, opener = Path(output), _opener(output)
                with opener(output, 'wt') as o:
                    o.write(s)
                return output
            else:
                return None
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
                    try:
                        score = float(line.split('=')[1].lstrip().split()[0])
                    except Exception as e:
                        logger.error(f'Failed to get docking score due to {e}:\n{line}')
                elif 'REMARK Name =' in line:
                    title = line.strip().split(' = ')[1]
        return score, title, ''.join(atom)

    def pdbqt(self, output=None, title=''):
        title = title or self.title
        if self.atom:
            s = [f'REMARK  Name = {title}' if title else '', '' if self.score is None else f'score {self.score}',
                 self.atom]
            s = '\n'.join([x for x in s if x])
            if output:
                with open(output, 'w') as o:
                    o.write(s)
                    return output
            else:
                return s
        else:
            return ''

    def sdf(self, output=None, title=''):
        _, pdbqt = tempfile.mkstemp(suffix='.pdbqt')
        _, _sdf = tempfile.mkstemp(suffix='.sdf')
        try:
            self.pdbqt(output=pdbqt, title=title)
            p = cmder.run(f'obabel {pdbqt} -o sdf -O {_sdf}', exit_on_error=False, fmt_cmd=False, log_cmd=False)
            if p.returncode:
                logger.error('Failed to write molecule to SDF file')
                return ''
            else:
                with open(_sdf) as f:
                    lines = (line for line in f
                             if not line.startswith('>  <REMARK>') and not line.startswith('  Name ='))
                    lines = ('\n' if 'OpenBabel' in line else line for line in lines)
                    if output:
                        with open(output, 'w') as o:
                            o.writelines(line for line in lines)
                        return output
                    else:
                        return ''.join(lines)
        finally:
            os.unlink(pdbqt)
            os.unlink(_sdf)

    def __str__(self):
        return self.s


def parse_sdf(sdf):
    path = str(sdf)
    opener = gzip.open if path.endswith('.sdf.gz') or path.endswith('.sdfgz') else open
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
    if path.endswith('.sdf.gz') or path.endswith('.sdfgz') or path.endswith('.sdf'):
        return parse_sdf(path)
    elif path.endswith('.dlg.gz') or path.endswith('.dlg'):
        return parse_dlg(path)
    else:
        raise TypeError(f'Unsupported file {path}')


def write(records, output=None, outdir=None):
    if output:
        opener = _opener(output)
        with opener(output, 'wt') as o:
            n = 0
            for record in records:
                if record:
                    o.write(record.sdf())
                    n += 1
            logger.debug(f'Successfully saved {n:,} items into {output}')
        return Path(output)
    elif outdir:
        outdir, outputs = Path(outdir), []
        outdir.mkdir(parents=True, exist_ok=True)
        for i, record in enumerate(records):
            if record:
                outputs.append(record.sdf(output=f'{outdir}/{record.title or i}.sdf'))
        logger.debug(f'Successfully saved {len(outputs):,} items into {outdir}')
    else:
        return ''.join(record for record in records if record)


def count_sdf(sdf):
    return sum(1 for m in parse_sdf(sdf) if m.mol)


def split_sdf(sdf, prefix, suffix='.sdf', files=0, records=0):
    names = []
    if files:
        n = count_sdf(sdf)
        num, remains = divmod(n, files)
        if remains:
            num += 1

        idx, end, items = 1, num, []
        for i, item in enumerate(parse_sdf(sdf), 1):
            items.append(item)
            if i == end:
                end += num
                name = write(items, output=f'{prefix}{idx}{suffix}')
                names.append(name)
                idx += 1
                items = []
        if items[:-1]:
            name = write(items[:-1], output=f'{prefix}{idx}{suffix}')
            names.append(name)
    elif records:
        idx, n, items = 0, records - 1, []
        for i, item in enumerate(parse_sdf(sdf), 0):
            items.append(item)
            idx, remains = divmod(i, records)

            if remains == n:
                name = write(items, output=f'{prefix}{idx + 1}{suffix}')
                names.append(name)
                items = []
        if items[:-1]:
            name = write(items[:-1], output=f'{prefix}{idx + 1}{suffix}')
            names.append(name)
    else:
        idx = 1
        name = write(parse_sdf(sdf), output=f'{prefix}{idx}{suffix}')
        names.append(name)
    return names


def sample_sdf(sdf, output, n=0, p=0.0, seed=None):
    logger.debug(f'Counting number of compounds in {sdf}')
    total = count_sdf(sdf)
    if n:
        logger.debug(f'Try to sample {n:,} compounds out of {total:,} compounds')
    elif p:
        n = int(total * p / 100)
        logger.debug(f'Try to sample {n:,} ({p}%) compounds out of {total:,} compounds')
    else:
        raise ValueError('Neither number of compounds (n) nor percentage of compounds (p) was specified, aborted!')

    if n <= total:
        logger.debug(f'Sampling {n:,} compounds from {sdf}')
    else:
        raise ValueError(f'No enough compounds ({total} < {n}) found in {sdf} to sample')

    logger.debug(f'Selecting random compounds')
    random.seed(seed)
    indices, ss = set(random.sample(range(total), n)), []
    for i, s in enumerate(parse_sdf(sdf)):
        if i in indices:
            ss.append(s)

    logger.debug(f'Saving {len(ss):,} compounds into {output}')
    opener = gzip.open if str(output).endswith('.sdf.gz') or str(output).endswith('.sdfgz') else open
    with opener(output, 'wt') as o:
        o.writelines(s.sdf() for s in ss if s.mol)
    logger.debug(f'Successfully saved {len(ss):,} compounds into {output}')

    # random.seed(seed)
    # indices = sorted(random.sample(range(total), n))
    # step = 1000
    # chunks = [indices[step * i: step * (i + 1)] for i in range(int(n / step) + 1)]
    #
    # ss, idx = [], 0
    # idx = 0
    # chunk, upper = chunks[idx], chunks[idx][-1]
    # for i, s in enumerate(parse_sdf(sdf)):
    #     if i < upper:
    #         if i in chunk:
    #             ss.append(s)
    #     elif i == upper:
    #         ss.append(s)
    #         idx += 1
    #         try:
    #             chunk, upper = chunks[idx], chunks[idx][-1]
    #         except IndexError:
    #             break
    #
    # opener = gzip.open if str(output).endswith('.sdf.gz') or str(output).endswith('.sdfgz') else open
    # with opener(output, 'wt') as o:
    #     o.writelines(s.sdf() for s in ss if s.mol)


def merge_sdf(sdfs, output, sort=None, max_score=0):
    items = []
    for sdf in sdfs:
        logger.debug(f'Parsing {sdf} ...')
        for s in parse_sdf(sdf):
            if s.score is not None and s.score < max_score:
                items.append(s)

    n = len(items)
    if sort:
        logger.debug(f'Sorting {n:,} docking poses on docking score {sort} ...')
        items = sorted(items, key=lambda x: x.score, reverse=False if sort == 'descending' else True)
        logger.debug(f'Sorting {n:,} docking poses on docking score {sort} complete.')

    write(items, output)
    logger.debug(f'Successfully saved {n:,} docking poses to {output}')
    return output


def batch_sdf(sdf, n, prefix):
    batches = [f'{prefix}{i + 1}.sdf' for i in range(n)]
    batches = [batch if Path(batch).exists() else '' for batch in batches]
    if all(batches):
        logger.debug('Ligand batches already exist, skip re-generate batches')
    else:
        logger.debug(f'Splitting ligands into {n:,} batches ...')
        batches = split_sdf(sdf, prefix, files=n)
        logger.debug(f'Successfully split ligands into {len(batches):,} batches')
    return batches


def clean_sdf(sdf, output=None):
    items = [sdf] if isinstance(sdf, SDF) else parse(sdf)
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


def dlg2sdf(dlg, sdf=None, title=''):
    s = ''.join(f'{item.sdf(title=title)}' for item in parse_dlg(dlg))
    if sdf is None:
        return s
    else:
        sdf = sdf or str(Path(dlg).with_suffix('.sdf'))
        with open(sdf, 'w') as o:
            o.write(s)
        return sdf


if __name__ == '__main__':
    pass
