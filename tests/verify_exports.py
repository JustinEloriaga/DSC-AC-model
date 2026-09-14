#!/usr/bin/env python3
"""Independently verify MATLAB data exports against the immutable source CSV.

Uses numpy/pandas, not the MATLAB preparation or correlation implementation.
Run from any directory; --data can point to the closure-masked sensitivity.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def verify(folder):
    summary = json.loads((folder / 'data_summary.json').read_text())
    source = Path(summary['source_file'])
    assert hashlib.sha256(source.read_bytes()).hexdigest() == summary['source_hash']
    tickers = summary['tickers']
    source_tickers = source.read_text().splitlines()[1].split(',')[1:]
    raw = pd.read_csv(source, skiprows=3, header=None, names=['Date'] + source_tickers)
    raw.index = pd.to_datetime(raw.pop('Date'), format='%d/%m/%Y')
    raw = raw.loc['2003-08-04':, tickers]
    assert len(tickers) == 14 and not any('USDCNH' in t for t in tickers)
    assert raw.index.is_unique and raw.index.is_monotonic_increasing
    Fridays = pd.date_range(raw.index.min(), raw.index.max(), freq='W-FRI')
    expected_levels = raw.reindex(raw.index.union(Fridays)).ffill().loc[Fridays]
    expected_returns = 100 * np.diff(np.log(expected_levels.to_numpy()), axis=0)
    levels = pd.read_csv(folder / 'weekly_levels.csv', index_col='Date')
    returns = pd.read_csv(folder / 'weekly_returns.csv', index_col='Date')
    assert list(levels.columns) == list(returns.columns) == tickers
    assert list(levels.index) == list(Fridays.strftime('%Y-%m-%d'))
    assert list(returns.index) == list(Fridays[1:].strftime('%Y-%m-%d'))
    assert returns.shape == (1204, 14) and levels.shape == (1205, 14)
    np.testing.assert_allclose(levels, expected_levels, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(returns, expected_returns, rtol=1e-9, atol=1e-10)

    quality = pd.read_csv(folder / 'weekly_quality.csv')
    for j, ticker in enumerate(tickers):
        observed = pd.Series(raw.index, index=raw.index).where(raw[ticker].notna())
        observed = observed.reindex(observed.index.union(Fridays)).ffill().loc[Fridays]
        q = quality.loc[quality.Ticker.eq(ticker)].sort_values('LevelDate')
        assert list(q.SourceDate) == list(observed.dt.strftime('%Y-%m-%d'))
        ages = (Fridays - pd.DatetimeIndex(observed)).days
        np.testing.assert_array_equal(q.AgeDays, ages)
    assert quality.Carried.sum() == 343
    assert quality.ClosureLevel.sum() == 1
    rq = pd.read_csv(folder / 'weekly_return_quality.csv')
    assert rq.ClosureReturn.sum() == 2
    x = returns.copy()
    for row in rq.loc[rq.ObservedForModel.eq(0)].itertuples():
        x.loc[row.Date, row.Ticker] = np.nan

    expected_pairs = [(i,j) for j in range(14) for i in range(j+1,14)]
    full = pd.read_csv(folder / 'full_sample_correlations.csv').sort_values('PairIndex')
    assert len(full) == 91
    max_corr_error = 0.0
    for (i,j), row in zip(expected_pairs,full.itertuples()):
        assert (row.TickerI,row.TickerJ) == (tickers[i],tickers[j])
        pair = x.iloc[:,[i,j]].dropna()
        rho = pair.corr().iloc[0,1]
        max_corr_error = max(max_corr_error,abs(rho-row.Correlation))
        assert row.N == len(pair)
        np.testing.assert_allclose(rho,row.Correlation,atol=1e-12)

    rolling = pd.read_csv(folder / 'rolling_correlations.csv')
    assert len(rolling) == 205114
    max_rolling_error = 0.0
    # Check every exported date/pair, not only selected charts.
    for (window,date), rows in rolling.groupby(['Window','Date'],sort=False):
        end = x.index.get_loc(date)+1
        block = x.iloc[end-window:end]
        assert len(block)==window
        c = block.corr().to_numpy()
        rows = rows.sort_values('PairIndex')
        assert list(rows.PairIndex)==list(range(1,92))
        for (i,j),row in zip(expected_pairs,rows.itertuples()):
            n = block.iloc[:,[i,j]].notna().all(axis=1).sum()
            assert row.N==n
            if n!=window:
                assert pd.isna(row.Correlation)
            else:
                err=abs(c[i,j]-row.Correlation)
                max_rolling_error=max(max_rolling_error,err)
                assert err<1e-12
    candidates=pd.read_csv(folder/'adf_candidates.csv')
    result={'status':'passed','returns':len(returns),'variables':len(tickers),
            'pairs':len(full),'rolling_rows':len(rolling),'source_hash':summary['source_hash'],
            'max_full_correlation_error':max_corr_error,
            'max_rolling_correlation_error':max_rolling_error,
            'adf_candidate_rows':len(candidates)}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=ROOT/'outputs/weekly_research/data')
    parser.add_argument('--json',type=Path,help='Optional destination for verification evidence')
    args=parser.parse_args()
    result=verify(args.data)
    message=json.dumps(result,indent=2)
    print(message)
    if args.json:
        args.json.parent.mkdir(parents=True,exist_ok=True)
        args.json.write_text(message+'\n')
