/* BOREXINO Reconstruction program
 * 
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: parameter_broker.hh,v 1.8 2006/12/05 13:37:20 razeto Exp $
 *
 * The parameter broker is a singleton containing the parameter
 * persistency: each parameter set by config files or by command 
 * line is stored here. This storage is exported to bx_named objects
 * which can fetch the parameter values (and so modules which are
 * bx_named). The parameter are indexed by bx_named name and parameter
 * name.
 * The interface is the follwoing: parameter can be fetched with get_parameter
 * which however throws an exception if the parameter is not present; a previus
 * check_parameter should be called to verify that the parameter is present.
 * The return value of get_parameter is a vtd object, see vtd.hh.
 *
 * Init parameter initialize a parameter, while set parameter checks that
 * the required parameter is existing before setting it; if not a warning is
 * generated.
 *
 */

#ifndef _PARAMETER_BROKER_H
#define _PARAMETER_BROKER_H

#include "vdt.hh"

#include "cmap.hh"
#include <vector>
#include <stdexcept>

class parameter_broker {
  public:
    parameter_broker (): parameter_matrix("parameter_matrix") {}

      // get is a bit long but a check is done for error in cmap. If error call internal method
    const vdt& get_parameter (const std::string& obj_name, const std::string& param_name) const { try { return parameter_matrix[obj_name][param_name]; } catch (std::runtime_error &er) { m_error_param_not_found (er, obj_name, param_name); } return unused; } 
    bool check_parameter (const std::string& obj_name, const std::string& param_name) const { return (parameter_matrix.check (obj_name) && parameter_matrix[obj_name].check(param_name)) ? true : false; }

    void init_parameter (const std::string& obj_name, const std::string& param_name, const std::string& value);
    void init_parameter (const std::string& obj_name, const std::string& param_name, const vdt& value);
    void set_parameter  (const std::string& obj_name, const std::string& param_name, const std::string& value);
    void set_parameter  (const std::string& obj_name, const std::string& param_name, const vdt& value);
    void set_parameter  (const std::string& obj_name_dot_param_name, const std::string& value);

    void set_configuration_name (const std::string name) { configuration_name = name; }
    const std::string get_configuration_name () const { return configuration_name; }

    void log_configuration () const;
  private:
    typedef std::cmap<std::string, vdt> parameter_namespace;
    std::cmap<std::string, parameter_namespace> parameter_matrix;
    std::string configuration_name;
    void m_error_param_not_found (const std::runtime_error &er, const std::string& obj_name, const std::string& param_name) const;
    vdt unused; // Really unused but usefull ;-P It is to return a dummy value on get_parameter in case of exception 
                // (and in that case it is not used).
};
#endif
/*
 * $Log: parameter_broker.hh,v $
 * Revision 1.8  2006/12/05 13:37:20  razeto
 * Added configuration name
 *
 * Revision 1.7  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.5  2004/11/24 13:04:33  razeto
 * Upgraded to the new cmap ctor with name
 *
 * Revision 1.4  2004/09/30 14:37:22  razeto
 * Added parameter_broker::log_configuration to log the current configuration
 *
 * Revision 1.3  2004/08/30 17:05:15  razeto
 * Upgraded to check for parameter existance error
 *
 * Revision 1.2  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.1  2004/06/10 12:38:19  razeto
 * Added bx_named and parameter_broker interface
 *
 */
