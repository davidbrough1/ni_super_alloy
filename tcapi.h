#include "tc_data_defs.h"

/*****************************
*   Appends the named database 
*/
void    tc_append_database(TC_STRING name);

/*****************************
*    Returns true if the checked license is available and valid, currently accepts: TC_DLL, TC_GUI and TC_TC4U
*/
TC_BOOL    tc_check_license(TC_STRING name,TC_STRING message,TC_STRING_LENGTH strlen_message);

/*****************************
*    Returns the status of “component_name” where status may be one of “ENTERED” or “SUSPENDED”
*/
TC_STRING    tc_component_status(TC_STRING component_name);

/*****************************
*    Computes the equilibrium in POLY-3 using the currently set conditions
*/
void    tc_compute_equilibrium();

/*****************************
*    Creates a new equilibrium in POLY-3 with number “equilibrium number”.
*/
void    tc_create_new_equilibrium(TC_INT equilibrium);

/*****************************
*    Returns the number of databases in the system and their names in “names_of_databases”
*/
TC_INT  tc_database(TC_STRING datan,TC_INT linelen);

/*****************************
*    Redefines the components in the system to the components in “component_name”.
*/
void    tc_define_components(TC_STRING component_name,TC_INT strlen,TC_INT number_of_components);

/*****************************
*    Returns the degrees of freedom in the system. This must be zero in order to perform an equilibrium calculation.
*/
TC_INT    tc_degrees_of_freedom();

/*****************************
*    
*/
void    tc_deinit();

/*****************************
*    Deletes the condition for the expression in “condition”.
*/
void    tc_delete_condition(TC_STRING condition);

/*****************************
*    Deletes a symbol in the system.
*/
void    tc_delete_symbol(TC_STRING symbol);

/*****************************
*    Returns the number of elements in the database and their names in “elements”
*/
TC_INT  tc_element(TC_STRING elements,TC_INT linelen);

/*****************************
*    Rejects “element name” in the currently selected database.
*/
void    tc_element_reject(TC_STRING element_name);

/*****************************
*    Selects “element name” in the currently selected database.
*/
void    tc_element_select(TC_STRING element_name);

/*****************************
*    Enters a parameter expression
*/
void    tc_enter_ges5_parameter(TC_STRING parameter,TC_STRING expression);

/*****************************
*    Enters a symbol in the system, the symbol type may be one of “CONSTANT”, “VARIABLE”, “FUNCTION” or “TABLE”,
*    “argument type” defines which of the following arguments will be used,
*    1 indicates the integer argument, 2 the double argument and 3 the string argument.
*/
void    tc_enter_symbol(TC_STRING symbol,TC_STRING type,TC_INT argument_type,TC_INT integer_argument,TC_FLOAT double_argument,TC_STRING string_argument);

/*****************************
*    Returns true if an error has been set,
*    returning the error number in “error number” and its corresponding message in “error message”
*/
TC_BOOL    tc_error(TC_INT *error_number,TC_STRING message,TC_INT strlen);

/*****************************
*    Sends a command to the GES5 module as defined in the argument “command”
*/
void    tc_ges5(TC_STRING command);

/*****************************
*    Executes the database command “GET_DATA”
*/
void    tc_get_data();

/*****************************
*    Returns Gm and the first derivatives with respect to site-fractions in “arr1”
*    and the second derivatives in “arr2” as GM.Y1.Y1, GM.Y1.Y2, GM.Y2.Y2, GM.Y1.Y3, GM.Y2.Y3 … GM.YN.YN
*/
void    tc_get_derivatives(TC_STRING phase_name,TC_FLOAT *arr1,TC_FLOAT *arr2);

/*****************************
*    Retrieves the expression of a parameter name
*/
void    tc_get_ges5_parameter(TC_STRING parameter,TC_STRING expression,TC_INT strlenExpression);

/*****************************
*    Retrieves the symbol or state variable value from the POLY-3 module.
*/
TC_FLOAT    tc_get_value(TC_STRING symbol);

/*****************************
*    Initializes the Thermo-Calc system, must be called prior to anything else.
*/
TC_INT    tc_init_root();

/*****************************
*    Initializes the Thermo-Calc system, must be called prior to anything else.
*       tmppath     Path to directory for log file
*       tcpath      Path to Thermo-Calc installation (Used to find databases)
*/
TC_INT    tc_init_root3(TC_STRING tmppath, TC_STRING tcpath);

/*****************************
*    Returns the number of components in the system and their names in “component_name”.
*/
TC_INT    tc_list_component(TC_STRING component_name,TC_INT strlen);

/*****************************
*    Returns the number of conditions set and their values in “conditions”
*/
TC_INT    tc_list_conditions(TC_STRING conditions,TC_INT strlen);

/*****************************
*    Returns the number of phases in the system and their names in “phase_name”.
*/
TC_INT    tc_list_phase(TC_STRING phase_name,TC_INT strlen);

/*****************************
*    Returns the number of species in the system and their names in “species_name”.
*/
TC_INT    tc_list_species(TC_STRING species_name,TC_INT strlen);

/*****************************
*    Returns the number of defined symbols in the system with their expression and
*    value in “symbols” and their corresponding type in “type of symbol”,
*    where the type may be one of
*    1=”CONSTANT”, 2=”VARIABLE”
*    3=”FUNCTION”
*    4=”TABLE”
*/
TC_INT    tc_list_symbols(TC_STRING symbols,TC_INT strlen,TC_INT *type);

/*****************************
*    
*/
TC_INT    tc_nr_of_constituents_in_phase(TC_STRING phase_name);

/*****************************
*    Opens the named database “name_of_database”
*/
void    tc_open_database(TC_STRING name);

/*****************************
*    Returns the number of phases in the system with the selected elements.
*    NOTE: the routine returns the number of all available phases.
*/
TC_INT    tc_phase(TC_STRING phases,TC_INT linelen);

/*****************************
*    Returns the number of sublattices in the phase (including phases with the status SUSPENDED),
*    the number of constituents on each sublattice in “constituents”,
*    the name of the selected species on each sublattice one after each other in “species names”
*    and the “number of sites” on each sublattice.
*/
TC_INT    tc_phase_all_constituents(TC_STRING phase_name,TC_INT *constituent_array,TC_STRING element_array,TC_INT strLenElem,TC_FLOAT *number_of_sites);

/*****************************
*    Returns the number of sublattices in the phase, the number of constituents on each sublattice in “constituents”,
*    the name of the selected species on each sublattice one after each other in “species names”
*    and the “number of sites” on each sublattice.
*/
TC_INT    tc_phase_constituents(TC_STRING phase_name,TC_INT *constituent_array,TC_STRING element_array,TC_INT strLenElem,TC_FLOAT *number_of_sites);

/*****************************
*    Rejects the phase in “phase_name”
*/
void    tc_phase_reject(TC_STRING phase_name);

/*****************************
*    Selects the phase in “phase_name”
*/
void    tc_phase_select(TC_STRING phase_name);

/*****************************
*    Returns the status of “phase_name” where status may be one of “FIXED”, “SUSPENDED” or “ENTERED”
*/
TC_STRING    tc_phase_status(TC_STRING phase_name);

/*****************************
*    Returns the number of sublattices in the phase, the number of constituents on each sublattice in “constituents”,
*    the name of the species on each sublattice one after each other in “species names”
*    and the number of sites in “number of sites”.
*/
TC_INT    tc_phase_structure(TC_STRING phase_name,TC_INT *constituent_array,TC_STRING species_array,TC_STRING_LENGTH strLenSpecies,TC_FLOAT *number_of_sites);

/*****************************
*    Sends a command to POLY-3 module as defined in the argument “command”
*/
void    tc_poly3(TC_STRING command);

/*****************************
*    
*/
void    tc_put_sitefractions(TC_STRING phase_name,TC_FLOAT *sfarr);

/*****************************
*    Loads the workspace from file “filename” in POLY-3.
*/
void    tc_read_poly3_file(TC_STRING filename);

/*****************************
*    Rejects the constituent “constituent” on sublattice “sublattiice” from phase “phase_name”.
*/
void    tc_reject_constituent(TC_STRING phase_name,TC_INT sublattice,TC_STRING constituent);

/*****************************
*    Resets the error if an error has been set
*/
void    tc_reset_error();

/*****************************
*    Restores the constituent “constituent” on sublattice “sublattiice” from phase “phase_name”.
*/
void    tc_restore_constituent(TC_STRING phase_name,TC_INT sublattice,TC_STRING constituent);

/*****************************
*    Stores/overwrites the current workspace in POLY-3 on the file “filename”.
*/
void    tc_save_poly3_file(TC_STRING filename);

/*****************************
*    Selects an equilibrium in POLY-3 with number “equilibrium number”.
*/
void    tc_select_equilibrium(TC_INT equilibrium);

/*****************************
*    Sets the status of “component_name” to “status” to one of “ENTERED” or “SUSPENDED”
*/
void    tc_set_component_status(TC_STRING component_name,TC_STRING status);

/*****************************
*    Sets a condition for the expression in “condition” to value in “value”.
*/
void    tc_set_condition(TC_STRING condition,TC_FLOAT value);

/*****************************
*    
*/
void    tc_set_license_code(TC_INT code);

/*****************************
*    Sets parameters for global minimization
*/
void    tc_set_minimization_option(TC_INT *global_flag,TC_INT *max_gridpoints,TC_INT *frequency,TC_INT *mesh_flag);

/*****************************
*    Sets the addition “addition” to the Gibbs
*/
void    tc_set_phase_addition(TC_STRING phase_name,TC_FLOAT addition);

/*****************************
*    Sets the status of “phase_name” to “status” to one of “FIXED”, “SUSPENDED”, DORMANT” or “ENTERED”.
*/
void    tc_set_phase_status(TC_STRING phase_name,TC_STRING status,TC_FLOAT value);

/*****************************
*    Sets a starting value for the “state variable” to “start value”.
*/
void    tc_set_start_value(TC_STRING state_variable,TC_FLOAT starting_value);

/*****************************
*    Returns the status of “species_name” where status may be one of “ENTERED” or “SUSPENDED”
*/
TC_STRING    tc_species_status(TC_STRING species_name);

/*****************************
*    Returns the version of Thermo-Calc in “version_name”
*/
void    tc_version(TC_STRING str,TC_INT str_len);
