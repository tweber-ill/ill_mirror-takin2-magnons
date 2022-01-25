/**
 * S(q,w) module for magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2022
 * @license GPLv2, see 'LICENSE' file
 */

// g++ -std=c++20 -I../ext -I../ext/takin -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -shared -fPIC -o magnonmod.so magnonmod.cpp ../ext/takin/tools/monteconvo/sqwbase.cpp ../ext/tlibs/log/log.cpp -llapacke

#include "magnonmod.h"

#include "libs/version.h"
#include "tlibs/string/string.h"
#include "tlibs/math/math.h"
#include "tlibs/phys/neutrons.h"


using t_real = typename MagnonMod::t_real;


// ----------------------------------------------------------------------------
// constructors

MagnonMod::MagnonMod()
{
	SqwBase::m_bOk = true;
}


MagnonMod::MagnonMod(const std::string& cfg_file) : MagnonMod()
{
	tl::log_info("Magnon module config file: \"", cfg_file, "\".");

	// load config file
	SqwBase::m_bOk = m_dyn.Load(cfg_file);
}


MagnonMod::~MagnonMod()
{
	m_dyn.Clear();
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// dispersion, spectral weight and structure factor

std::tuple<std::vector<t_real>, std::vector<t_real>>
	MagnonMod::disp(t_real h, t_real k, t_real l) const
{
	// calculate dispersion relation
	std::vector<t_real> energies;
	std::vector<t_real> weights;

	auto modes = m_dyn.GetEnergies(h, k, l, false);

	for(const auto& mode : modes)
	{
		energies.push_back(mode.E);
		weights.push_back(mode.weight);
	}

	return std::make_tuple(energies, weights);
}


t_real MagnonMod::operator()(t_real h, t_real k, t_real l, t_real E) const
{
	std::vector<t_real> Es, Ws;
	std::tie(Es, Ws) = disp(h, k, l);

	t_real incoh = 0.;
	if(!tl::float_equal(m_incoh_amp, t_real(0)))
		incoh = tl::gauss_model(E, t_real(0),
			m_incoh_sigma, m_incoh_amp, t_real(0));

	t_real S = 0.;
	for(std::size_t iE=0; iE<Es.size(); ++iE)
	{
		if(!tl::float_equal(Ws[iE], t_real(0)))
			S += tl::gauss_model(E, Es[iE], m_sigma, Ws[iE], t_real(0));
	}

	return m_S0*S /** tl::bose_cutoff(E, m_T, cutoff)*/ + incoh;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// get & set variables

std::vector<MagnonMod::t_var> MagnonMod::GetVars() const
{
	std::vector<t_var> vars;

	// get bragg peak
	const auto& _G = m_dyn.GetBraggPeak();
	std::vector<t_real> G;
	if(_G.size() == 3)
	{
		G = std::vector<t_real>{{
			_G[0].real(),
			_G[1].real(),
			_G[2].real()
		}};
	}
	else
	{
		G = std::vector<t_real>{{ 1., 0., 0. }};
	}

	// get external field
	const tl2_mag::ExternalField& field = m_dyn.GetExternalField();
	std::vector<t_real> B;
	if(field.dir.size() == 3)
	{
		B = std::vector<t_real>{{
			field.dir[0].real(),
			field.dir[1].real(),
			field.dir[2].real()
		}};
	}
	else
	{
		B = std::vector<t_real>{{ 0., 0., 1. }};
	}

	vars.push_back(SqwBase::t_var{
		"sigma", "real", tl::var_to_str(m_sigma)});
	vars.push_back(SqwBase::t_var{
		"inc_amp", "real", tl::var_to_str(m_incoh_amp)});
	vars.push_back(SqwBase::t_var{
		"inc_sigma", "real", tl::var_to_str(m_incoh_sigma)});
	vars.push_back(SqwBase::t_var{
		"S0", "real", tl::var_to_str(m_S0)});
	vars.push_back(SqwBase::t_var{
		"T", "real", tl::var_to_str(m_dyn.GetTemperature())});
	vars.push_back(SqwBase::t_var{
		"cutoff", "real", tl::var_to_str(m_dyn.GetBoseCutoffEnergy())});
	vars.push_back(SqwBase::t_var{
		"G", "vector", vec_to_str(G)});
	vars.push_back(SqwBase::t_var{
		"B_dir", "vector", vec_to_str(B)});
	vars.push_back(SqwBase::t_var{
		"B_mag", "real", tl::var_to_str(field.mag)});
	vars.push_back(SqwBase::t_var{
		"B_align_spins", "real", tl::var_to_str((int)field.align_spins)});

	return vars;
}


void MagnonMod::SetVars(const std::vector<MagnonMod::t_var>& vars)
{
	if(!vars.size()) return;

	for(const SqwBase::t_var& var : vars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "sigma")
			m_sigma = tl::str_to_var<t_real>(strVal);
		else if(strVar == "inc_amp")
			m_incoh_amp = tl::str_to_var<decltype(m_incoh_amp)>(strVal);
		else if(strVar == "inc_sigma")
			m_incoh_sigma = tl::str_to_var<decltype(m_incoh_sigma)>(strVal);
		else if(strVar == "S0")
			m_S0 = tl::str_to_var<decltype(m_S0)>(strVal);
		else if(strVar == "T")
			m_dyn.SetTemperature(tl::str_to_var<t_real>(strVal));
		else if(strVar == "cutoff")
			m_dyn.SetBoseCutoffEnergy(tl::str_to_var<t_real>(strVal));
		else if(strVar == "G")
		{
			std::vector<t_real> G = str_to_vec<std::vector<t_real>>(strVal);
			if(G.size() == 3)
			{
				m_dyn.SetBraggPeak(G[0], G[1], G[2]);
				m_dyn.CalcSpinRotation();
			}
			else
			{
				tl::log_err("Invalid Bragg peak.");
			}
		}
		else if(strVar == "B_dir")
		{
			std::vector<t_real> dir = str_to_vec<std::vector<t_real>>(strVal);
			if(dir.size() == 3)
			{
				tl2_mag::ExternalField field = m_dyn.GetExternalField();
				field.dir = tl2::create<tl2_mag::t_vec>({
					dir[0], dir[1], dir[2]});
				m_dyn.SetExternalField(field);
				m_dyn.CalcSpinRotation();
			}
			else
			{
				tl::log_err("Invalid field direction.");
			}
		}
		else if(strVar == "B_mag")
		{
			tl2_mag::ExternalField field = m_dyn.GetExternalField();
			field.mag = tl::str_to_var<decltype(m_S0)>(strVal);
			m_dyn.SetExternalField(field);
			m_dyn.CalcSpinRotation();
		}
		else if(strVar == "B_align_spins")
		{
			tl2_mag::ExternalField field = m_dyn.GetExternalField();
			field.align_spins = (tl::str_to_var<int>(strVal) != 0);
			m_dyn.SetExternalField(field);
			m_dyn.CalcSpinRotation();
		}
	}
}


bool MagnonMod::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	return SqwBase::SetVarIfAvail(strKey, strNewVal);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// copy

SqwBase* MagnonMod::shallow_copy() const
{
	MagnonMod *mod = new MagnonMod();

	mod->m_sigma = this->m_sigma;
	mod->m_incoh_amp = this->m_incoh_amp;
	mod->m_incoh_sigma = this->m_incoh_sigma;
	mod->m_S0 = this->m_S0;
	mod->m_dyn = this->m_dyn;

	return mod;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// SO interface

#include <boost/dll/alias.hpp>


std::tuple<std::string, std::string, std::string> sqw_info()
{
	tl::log_info("In ", __func__, ".");
	return std::make_tuple(TAKIN_VER, "magnonmod", "Magnon Dynamics");
}


std::shared_ptr<SqwBase> sqw_construct(const std::string& cfg_file)
{
	tl::log_info("In ", __func__, ".");
	return std::make_shared<MagnonMod>(cfg_file);
}


// exports from so file
BOOST_DLL_ALIAS(sqw_info, takin_sqw_info);
BOOST_DLL_ALIAS(sqw_construct, takin_sqw);
// ----------------------------------------------------------------------------
