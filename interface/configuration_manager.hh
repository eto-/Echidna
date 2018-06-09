/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: configuration_manager.hh,v 1.5 2006/12/05 12:17:18 razeto Exp $
 *
 * The configuration_manager is the object which collects every
 * configuration bit from files and command line (parameters). When the
 * `process of collecting is finished it sends the module list to the module_factory
 * and the parameters of the current configuration to the parameter_broker.
 * 
 * For parsing the configuration files a set of handler classes are present, but
 * to avoid external usage they are in the class private namespace.
 *
 */
#ifndef _CONFIGURATION_MANAGER_H
#define _CONFIGURATION_MANAGER_H
#include "messenger.hh"
#include "cmap.hh"
#include "parameter_broker.hh"

#include <string>
#include <vector>
#include <iostream>

class configuration_manager: public parameter_broker {
  public:
    configuration_manager () {}
    
    void add_main_config_file (const std::string& filename);
    void add_user_config_file (const std::string& filename);
    void add_comman_line_parameter (const std::string& name_dot_parameter, const std::string& value);

    std::vector<std::string> get_module_name_list () const;
    void upload_configuration (const std::string& configuration_name);

    class configuration_parameter;
    class configuration_namespace;
    void add_item (const configuration_namespace& item);

    const configuration_namespace& get_configuration (const std::string& configuration_name);
  private:
    std::vector<configuration_namespace> module_list, named_list, config_list, user_file_params, command_line_params;

  public:

      // The configuration_namespace is the collection of related parameters: it can represent a module
      // description, a configuration or the set of command line parameters.
      // It contains a vector of configuration_parameter.
      // The upload_parameter_broker method for each element of the configuration_parameter vector upload the
      // parameter broker with the element values: if this is module or named the element is considered part
      // of a module/named initialization and so the parameter_broker::init() is called for uploading.
      // Alternativelly config, user_file and command_line_parameters should be a list of assigmenent statement,
      // and the configuration_parameter name should be module_name.parameter_name (or named.parameter). This case 
      // parameter_broker::set () is called.
    class configuration_namespace {
      public:
        enum item_type {
          module,
          named,
          config,
          user_file,
          command_line_parameters,
        };
        configuration_namespace (item_type type, const std::string& ns_name = ""): i_type(type), s_ns_name(ns_name) {}

        item_type get_type () const { return i_type; }
	const std::string& get_namespace_name () const { return s_ns_name; }

        void add_parameter (const configuration_parameter& param) { param_v.push_back (param); }
        void upload_parameter_broker (parameter_broker *broker) const;
	void operator+= (const configuration_namespace& right) { param_v.insert (param_v.end (), right.param_v.begin (), right.param_v.end ()); }
	
      private:
        item_type i_type;
        std::string s_ns_name;
        std::vector<configuration_parameter> param_v;
    };

      // A configuration_parameter is just the pair of name/value of a parameter; the name can be 
      // module_name.parameter_name or just parameter_name dependig of the semantic of the enclosing 
      // configuration_namespace.
    class configuration_parameter {
      public:
        configuration_parameter (const std::string& name, const std::string& value): s_name (name), s_value (value) {}
        const std::string& get_name () const { return s_name; }
        const std::string& get_value () const { return s_value; }
      private:
        std::string s_name, s_value;
    };
  private:
      // A small functor just to search configuration_namespace by names
    struct check_ns_name {
      check_ns_name (const std::string& ns_name): s_ns_name(ns_name) {}
      bool operator() (const configuration_namespace& item) { return item.get_namespace_name () == s_ns_name; }
      const std::string s_ns_name;
    };

  friend std::ostream& operator<< (std::ostream& out, const configuration_parameter& param);
};

inline std::ostream& operator<< (std::ostream& out, const configuration_manager::configuration_parameter& param) { 
  return out << "\"" << param.get_name () << " = " << param.get_value () << "\"";
}
#endif
/*
 * $Log: configuration_manager.hh,v $
 * Revision 1.5  2006/12/05 12:17:18  razeto
 * Add <include configuration> inside an other configuration
 *
 * Revision 1.4  2005/08/02 14:01:49  razeto
 * Removed some useless include
 *
 * Revision 1.3  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 */
