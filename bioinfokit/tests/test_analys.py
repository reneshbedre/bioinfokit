from bioinfokit.analys import norm, get_data, genfam, stat
import pytest
import numpy as np
import pandas as pd
from unittest import TestCase


class TestNormalization(TestCase):
    def test_cpm(self):
        df = get_data('sc_exp').data
        df = df.drop(['length'], axis=1)
        df = df.set_index('gene')
        nm = norm()
        nm.cpm(df=df)
        np.testing.assert_array_equal(round(nm.cpm_norm.iloc[0], 2).to_numpy(),
                                      np.asarray([100.695004, 101.731189, 74.721094, 92.633828, 74.270713, 95.314714])
                                      .round(2))

    def test_rpkm(self):
        df = get_data('sc_exp').data
        df = df.set_index('gene')
        nm = norm()
        nm.rpkm(df=df, gl='length')
        np.testing.assert_array_equal(round(nm.rpkm_norm.iloc[0], 2).to_numpy(),
                                      np.asarray([50.804745, 51.327542, 37.699846, 46.737552, 37.472610, 48.090169])
                                      .round(2))

    def test_tpm(self):
        df = get_data('sc_exp').data
        df = df.set_index('gene')
        nm = norm()
        nm.tpm(df=df, gl='length')
        np.testing.assert_array_equal(round(nm.tpm_norm.iloc[0], 2).to_numpy(),
                                      np.asarray([99.730156, 97.641941, 72.361658, 89.606265, 69.447237, 90.643338])
                                      .round(2))

    def test_genfam(self):
        id_data = pd.read_csv('https://reneshbedre.github.io/assets/posts/genfam/grai_id.txt', header=None)
        res = genfam()
        res.fam_enrich(id_file=id_data, species='grai', id_type=1)
        np.testing.assert_array_equal(res.df_enrich.iloc[0][1], 'MYB')

    '''
    def test_tukeyhsd(self):
        d = pd.read_csv("https://reneshbedre.github.io/assets/posts/anova/twowayanova.txt", sep="\t")
        d_melt = pd.melt(d, id_vars=['Genotype'], value_vars=['1_year', '2_year', '3_year'])
        d_melt.columns = ['Genotype', 'years', 'value']
        res = stat()
        res.tukey_hsd(df=d_melt, res_var='value', xfac_var=['Genotype', 'years'],
                      anova_model='value ~ C(Genotype) + C(years) + C(Genotype):C(years)')
        print(res.tukey_summary.head())
        np.testing.assert_array_equal(round(res.tukey_summary[(res.tukey_summary['group1'] == ('A', '1_year')) & (res.tukey_summary['group2'] == ('A', '2_year'))].iloc[0]['Diff'], 2), 2.38)
        np.testing.assert_array_equal(round(res.tukey_summary[(res.tukey_summary['group1'] == ('A', '1_year')) & (res.tukey_summary['group2'] == ('A', '2_year'))].iloc[0]['q-value'], 2), 6.89)
    '''