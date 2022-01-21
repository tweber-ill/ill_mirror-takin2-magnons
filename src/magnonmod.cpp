/**
 * S(q,w) module for magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2022
 * @license GPLv2, see 'LICENSE' file
 */

// g++ -std=c++20 -I -I../ext -I../ext/takin -shared -fPIC -o magnonmod.so magnonmod.cpp ../ext/takin/tools/monteconvo/sqwbase.cpp ../ext/tlibs/log/log.cpp

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
	SqwBase::m_bOk = true;
}


MagnonMod::~MagnonMod()
{
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// dispersion, spectral weight and structure factor

std::tuple<std::vector<t_real>, std::vector<t_real>>
	MagnonMod::disp(t_real h, t_real k, t_real l) const
{
	// TODO: calculate dispersion relation

	return std::make_tuple(
		std::vector<t_real>({0.}),  // energies
		std::vector<t_real>({0.})); // weights
}


t_real MagnonMod::operator()(t_real h, t_real k, t_real l, t_real E) const
{
	t_real cutoff = t_real(0.02);

	std::vector<t_real> vecE, vecW;
	std::tie(vecE, vecW) = disp(h, k, l);

	t_real incoh=0, S_p=0, S_m=0;
	if(!tl::float_equal(m_incoh_amp, t_real(0)))
		incoh = tl::gauss_model(E, t_real(0), m_incoh_sigma, m_incoh_amp, t_real(0));

	t_real dS = 0;
	for(std::size_t iE=0; iE<vecE.size(); ++iE)
	{
		if(!tl::float_equal(vecW[iE], t_real(0)))
			dS += tl::gauss_model(E, vecE[iE], m_sigma, vecW[iE], t_real(0));
	}

	return m_S0*dS * tl::bose_cutoff(E, m_T, cutoff) + incoh;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// get & set variables

std::vector<MagnonMod::t_var> MagnonMod::GetVars() const
{
	std::vector<t_var> vars;

	vars.push_back(SqwBase::t_var{"sigma", "real", tl::var_to_str(m_sigma)});
	vars.push_back(SqwBase::t_var{"inc_amp", "real", tl::var_to_str(m_incoh_amp)});
	vars.push_back(SqwBase::t_var{"inc_sigma", "real", tl::var_to_str(m_incoh_sigma)});
	vars.push_back(SqwBase::t_var{"S0", "real", tl::var_to_str(m_S0)});
	vars.push_back(SqwBase::t_var{"T", "real", tl::var_to_str(m_T)});

	return vars;
}


void MagnonMod::SetVars(const std::vector<MagnonMod::t_var>& vars)
{
	if(!vars.size()) return;

	for(const SqwBase::t_var& var : vars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "sigma") m_sigma = tl::str_to_var<t_real>(strVal);
		else if(strVar == "inc_amp") m_incoh_amp = tl::str_to_var<decltype(m_incoh_amp)>(strVal);
		else if(strVar == "inc_sigma") m_incoh_sigma = tl::str_to_var<decltype(m_incoh_sigma)>(strVal);
		else if(strVar == "S0") m_S0 = tl::str_to_var<decltype(m_S0)>(strVal);
		else if(strVar == "T") m_T = tl::str_to_var<decltype(m_T)>(strVal);
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
	mod->m_T = this->m_T;

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
