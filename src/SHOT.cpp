/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE
   This software is licensed under the Eclipse Public License 2.0.
   Please see the README and LICENSE files for more information.
*/

#include "Environment.h"
#include "Solver.h"
#include "Report.h"
#include "Utilities.h"
#include "Output.h"
#include "Settings.h"
#include "Problem.h"

#include "argh.h"

#ifdef HAS_STD_FILESYSTEM
#include <filesystem>
namespace fs = std;
#endif

#ifdef HAS_STD_EXPERIMENTAL_FILESYSTEM
#include <experimental/filesystem>
namespace fs = std::experimental;
#endif

#include <iostream>
#include <string>

using namespace SHOT;

int main(int argc, char* argv[])
{
    Solver solver;
    auto env = solver.getEnvironment();
    bool useASL = false;

    argh::parser cmdl;
    cmdl.add_params({ "--opt", "--osol" });
    cmdl.add_params({ "--osrl", "--trc", "--log" });
    cmdl.add_params({ "--sol" });
    cmdl.add_params({ "--docs" });
    cmdl.add_params({ "--debug" });

    cmdl.parse(argc, argv);

    std::string filename;
    fs::filesystem::path resultFile, optionsFile, traceFile, logFile, solFile;

    // Read or create the file for the log
    if(cmdl("--log") >> filename) // Have specified a log-file
    {
        logFile = fs::filesystem::current_path() / fs::filesystem::path(filename);
        solver.setLogFile(logFile.string());
    }
    else
    {
        logFile = fs::filesystem::current_path() / fs::filesystem::path("SHOT.log");
        solver.setLogFile(logFile.string());
    }

    env->report->outputSolverHeader();

    if(cmdl["--help"])
    {
        env->output->outputInfo("╶──────────────────────────────────────────────────────────────────────────────────"
                                "───────────────────────────────────╴\r\n");

        env->output->outputCritical("  Usage: SHOT PROBLEMFILE [ARGUMENTS] [OPTIONS]");
        env->output->outputCritical("");
        env->output->outputCritical("  SHOT has been compiled with support for the following problem formats ");

#ifdef HAS_AMPL
        env->output->outputCritical("   AMPL (.nl) ");
#endif
#ifdef HAS_GAMS
        env->output->outputCritical("   GAMS (.gms) ");
#endif
        env->output->outputCritical("   OSiL (.osil or .xml) ");
        env->output->outputCritical("");
        env->output->outputCritical("  The following command line arguments can also be used:");
        env->output->outputCritical("");
#ifdef HAS_AMPL
        env->output->outputCritical("   --AMPL                   Activates ASL support. Only to be used with nl-files");
#endif
        env->output->outputCritical("   --debug [DIRECTORY]      Saves debug information in the specified directory");
        env->output->outputCritical(
            "                            If DIRECTORY is empty 'debug/<problem_name>/' is used");
        env->output->outputCritical("   --log FILE               Sets the filename for the log file");
        env->output->outputCritical("   --opt [FILE]             Reads in options from FILE in GAMS format");
        env->output->outputCritical("                            If FILE is empty, a new options file will be created");
        env->output->outputCritical("   --osol [FILE]            Reads in options from FILE in OSoL format");
        env->output->outputCritical("                            If FILE is empty, a new options file will be created");
        env->output->outputCritical("   --osrl FILE              Sets the filename for the OSrL result file");
        env->output->outputCritical(
            "   --trc [FILE]             Prints a trace file to <problemname>.trc or specified filename");
        env->output->outputCritical("");
        env->output->outputCritical("");
        env->output->outputCritical("  It is possible to specify options directly using the the command line:");
        env->output->outputCritical("");
        env->output->outputCritical("   OPTIONNAME=VALUE OPTIONNAME=VALUE ...");
        env->output->outputCritical("");
        env->output->outputCritical("  See the documentation for valid options. For example:");
        env->output->outputCritical("");
        env->output->outputCritical("   Termination.TimeLimit=100.0   Sets the time limit to 100 seconds");
        env->output->outputCritical("");
        env->output->outputCritical("  Can also use the following shorthand options:");
        env->output->outputCritical("");
        env->output->outputCritical("   --convex                 Assumes the problem is convex");
#ifdef HAS_CBC
        env->output->outputCritical("   --mip=cbc                Sets the MIP solver to Cbc");
#endif
#ifdef HAS_CPLEX
        env->output->outputCritical("   --mip=cplex              Sets the MIP solver to Cplex");
#endif
#ifdef HAS_GUROBI
        env->output->outputCritical("   --mip=gurobi             Sets the MIP solver to Gurobi");
#endif
#ifdef HAS_GAMS
        env->output->outputCritical("   --nlp={ipopt, gams}      Sets the NLP solver to use");
#endif
#ifdef HAS_CPLEX
        env->output->outputCritical("   --tree={single, multi}   Activates single- or multi-tree strategy");
#elif HAS_GUROBI

        env->output->outputCritical("   --tree={single, multi}   Activates single- or multi-tree strategy");
#endif
        env->output->outputCritical("   --threads=VALUE          Sets the maximum number of threads to use");
        env->output->outputCritical("   --absgap=VALUE           Sets the absolute objective gap tolerance");
        env->output->outputCritical("   --relgap=VALUE           Sets the relative objective gap tolerance");
        env->output->outputCritical("   --timelimit=VALUE        Sets the time-limit in seconds");
        env->output->outputCritical("");

        return (0);
    }

    // Generate a markup file with the options
    if(cmdl["--docs"])
    {
        std::string markup = env->settings->getSettingsAsMarkup();

        auto filepath = fs::filesystem::current_path() / fs::filesystem::path("options.md");
        if(!Utilities::writeStringToFile(filepath.string(), markup))
            env->output->outputCritical(" Error when writing markup file: " + filepath.string());
    }

    // Read or create options file

    bool defaultOptionsGenerated = false;

    if(cmdl("--opt") >> filename) // Have specified a opt-file
    {
        auto filepath = fs::filesystem::current_path() / fs::filesystem::path(filename);

        if(fs::filesystem::exists(filepath))
        {
            optionsFile = filepath;
        }
        else
        {
            env->output->outputCritical("  Options file not found: " + filepath.string());
            return 0;
        }
    }
    else if(cmdl["--opt"]) // Create a new opt-file
    {
        auto filepath = fs::filesystem::current_path() / fs::filesystem::path("options.opt");

        if(fs::filesystem::exists(filepath))
        {
            optionsFile = filepath;
        }
        else
        {
            // Create option file
            if(!Utilities::writeStringToFile(filepath.string(), solver.getOptions()))
            {
                env->output->outputCritical("  Error when writing options file: " + filepath.string());
                return 0;
            }

            defaultOptionsGenerated = true;
            env->output->outputInfo("  Default options file written to: " + filepath.string());
        }
    }
    else if(cmdl("--osol") >> filename) // Have specified a OSoL-file
    {
        auto filepath = fs::filesystem::current_path() / fs::filesystem::path(filename);

        if(fs::filesystem::exists(filepath))
        {
            optionsFile = filepath;
            env->output->outputInfo("  Default options file written to: " + filepath.string());
        }
        else
        {
            env->output->outputCritical(" Options file not found: " + filepath.string());
            return 0;
        }
    }
    else if(cmdl["--osol"]) // Create a new OSoL-file
    {
        auto filepath = fs::filesystem::current_path() / fs::filesystem::path("options.xml");

        if(fs::filesystem::exists(filepath))
        {
            optionsFile = filepath;
        }
        else
        {
            // Create OSoL-file
            if(!Utilities::writeStringToFile(filepath.string(), solver.getOptionsOSoL()))
            {
                env->output->outputCritical("  Error when writing OSoL file: " + filepath.string());
                return 0;
            }
        }

        defaultOptionsGenerated = true;
    }

    if(!defaultOptionsGenerated)
    {
        if(!optionsFile.empty() && !solver.setOptionsFromFile(optionsFile.string()))
        {
            env->output->outputCritical("  Cannot set options from file: " + optionsFile.string());
            return (0);
        }

        env->output->setLogLevels(static_cast<E_LogLevel>(env->settings->getSetting<int>("Console.LogLevel", "Output")),
            static_cast<E_LogLevel>(env->settings->getSetting<int>("File.LogLevel", "Output")));
    }

    // Reads options specified in the command line arguments
#if HAS_AMPL
    if(cmdl["--AMPL"])
    {
        useASL = true;

        // We always want to write to where the problem is when called by ASL
        solver.updateSetting("OutputDirectory", "Output", static_cast<int>(ES_OutputDirectory::Problem));
    }
#endif

    if(cmdl["--convex"])
        solver.updateSetting("AssumeConvex", "Convexity", true);

    if(cmdl["--debug"])
        solver.updateSetting("Debug.Enable", "Output", true);

    std::string debugPath;
    if(cmdl("--debug") >> debugPath)
    {
        solver.updateSetting("Debug.Enable", "Output", true);
        solver.updateSetting("Debug.Path", "Output", debugPath);
    }

    std::string argValue;

    if(cmdl("--mip") >> argValue)
    {
#ifdef HAS_CBC
        if(argValue == "cbc")
            solver.updateSetting("MIP.Solver", "Dual", static_cast<int>(ES_MIPSolver::Cbc));
#endif
#ifdef HAS_CPLEX
        if(argValue == "cplex")
            solver.updateSetting("MIP.Solver", "Dual", static_cast<int>(ES_MIPSolver::Cplex));
#endif
#ifdef HAS_GUROBI
        if(argValue == "gurobi")
            solver.updateSetting("MIP.Solver", "Dual", static_cast<int>(ES_MIPSolver::Gurobi));
#endif
    }

    if(cmdl("--nlp") >> argValue)
    {
#ifdef HAS_GAMS
        if(argValue == "gams")
            solver.updateSetting("FixedInteger.Solver", "Primal", static_cast<int>(ES_PrimalNLPSolver::GAMS));
#endif
#ifdef HAS_IPOPT
        if(argValue == "ipopt")
            solver.updateSetting("FixedInteger.Solver", "Primal", static_cast<int>(ES_PrimalNLPSolver::Ipopt));
#endif
    }

    if(cmdl("--tree") >> argValue)
    {
#ifdef HAS_CPLEX
        if(argValue == "single")
            solver.updateSetting("TreeStrategy", "Dual", static_cast<int>(ES_TreeStrategy::SingleTree));
        else if(argValue == "multi")
            solver.updateSetting("TreeStrategy", "Dual", static_cast<int>(ES_TreeStrategy::MultiTree));
#endif
#ifdef HAS_GUROBI
        if(argValue == "single")
            solver.updateSetting("TreeStrategy", "Dual", static_cast<int>(ES_TreeStrategy::SingleTree));
        else if(argValue == "multi")
            solver.updateSetting("TreeStrategy", "Dual", static_cast<int>(ES_TreeStrategy::MultiTree));
#endif
    }

    if(cmdl("--threads") >> argValue)
    {
        try
        {
            solver.updateSetting("MIP.NumberOfThreads", "Dual", std::stoi(argValue));
        }
        catch(const std::exception& e)
        {
            env->output->outputCritical("  Cannot read value for parameter 'threads'");
        }
    }

    if(cmdl("--absgap") >> argValue)
    {
        try
        {
            solver.updateSetting("ObjectiveGap.Absolute", "Termination", std::stod(argValue));
        }
        catch(const std::exception& e)
        {
            env->output->outputCritical("  Cannot read value for parameter 'absgap'");
        }
    }

    if(cmdl("--relgap") >> argValue)
    {
        try
        {
            solver.updateSetting("ObjectiveGap.Relative", "Termination", std::stod(argValue));
        }
        catch(const std::exception& e)
        {
            env->output->outputCritical("  Cannot read value for parameter 'relgap'");
        }
    }

    if(cmdl("--timelimit") >> argValue)
    {
        try
        {
            solver.updateSetting("TimeLimit", "Termination", std::stod(argValue));
        }
        catch(const std::exception& e)
        {
            env->output->outputCritical("  Cannot read value for parameter 'timelimit'");
        }
    }

    for(auto& ARG : cmdl.pos_args())
    {
        int dotLocation = ARG.find('.');

        if(dotLocation == std::string::npos)
            continue;

        int equalLocation = ARG.find('=');

        if(equalLocation == std::string::npos)
            continue;

        if(equalLocation <= dotLocation)
            continue;

        auto category = ARG.substr(0, dotLocation);
        auto name = ARG.substr(dotLocation + 1, equalLocation - dotLocation - 1);
        auto value = ARG.substr(equalLocation + 1, ARG.size());

        bool found = false;

        for(auto& S : solver.getSettingIdentifiers(E_SettingType::String))
        {
            if(ARG.find(S) == 0)
            {
                solver.updateSetting(name, category, value);
                found = true;
                break;
            }
        }

        if(found)
            continue;

        for(auto& S : solver.getSettingIdentifiers(E_SettingType::Boolean))
        {
            if(ARG.find(S) == 0)
            {
                if(value == "true" || value == "false")
                    solver.updateSetting(name, category, (value == "true" ? true : false));
                else
                    env->output->outputCritical("  Cannot read boolean option in " + ARG);

                found = true;
                break;
            }
        }

        if(found)
            continue;

        for(auto& S : solver.getSettingIdentifiers(E_SettingType::Integer))
        {
            if(ARG.find(S) == 0)
            {
                try
                {
                    solver.updateSetting(name, category, std::stoi(value));
                }
                catch(const std::exception& e)
                {
                    env->output->outputCritical("  Cannot read integer option in " + ARG);
                }

                found = true;
                break;
            }
        }

        if(found)
            continue;

        for(auto& S : solver.getSettingIdentifiers(E_SettingType::Enum))
        {
            if(ARG.find(S) == 0)
            {
                try
                {
                    solver.updateSetting(name, category, std::stoi(value));
                }
                catch(const std::exception& e)
                {
                    env->output->outputCritical("  Cannot read integer option in " + ARG);
                }

                found = true;
                break;
            }
        }

        if(found)
            continue;

        for(auto& S : solver.getSettingIdentifiers(E_SettingType::Double))
        {
            if(ARG.find(S) == 0)
            {
                try
                {
                    solver.updateSetting(name, category, std::stod(value));
                }
                catch(const std::exception& e)
                {
                    env->output->outputCritical("  Cannot read numeric option in " + ARG);
                }

                found = true;
                break;
            }
        }
    }

    // Read problem file

    if(!cmdl(1) || !(cmdl(1) >> filename))
    {
        env->output->outputCritical("  No problem file specified.\r\n");
        env->output->outputCritical("  Try 'SHOT --help' for more information.");
        return (0);
    }

    if(!fs::filesystem::exists(filename))
    {
        if(useASL && fs::filesystem::exists(filename + ".nl"))
        {
            filename += ".nl";
        }
        else
        {
            env->output->outputCritical("   Problem file " + filename + " not found!");
            return (0);
        }
    }

    if(!solver.setProblem(filename))
    {
        env->output->outputCritical("   Error when reading problem file.");
        return (0);
    }

    // Check if we want to use the ASL calling format
    if(useASL && !((ES_SourceFormat)env->settings->getSetting<int>("SourceFormat", "Input") == ES_SourceFormat::NL))
    {
        env->output->outputCritical("  Error: Can only use parameter AMPL if the problem is a AMPL (.nl) file.");
        return (0);
    }

    // Define result file locations

    std::string osrlFilename;
    if(cmdl("--osrl") >> osrlFilename) // Have specified a OSrL-file location
    {
        resultFile = fs::filesystem::path(env->settings->getSetting<std::string>("ResultPath", "Output"))
            / fs::filesystem::path(osrlFilename);
    }

    std::string trcFilename;
    if(cmdl("--trc") >> trcFilename) // Have specified a trace-file location
    {
        traceFile = fs::filesystem::path(env->settings->getSetting<std::string>("ResultPath", "Output"))
            / fs::filesystem::path(trcFilename);
    }

    std::string solFilename;
    if(cmdl("--sol") >> solFilename) // Have specified an sol-file location
    {
        solFile = fs::filesystem::path(env->settings->getSetting<std::string>("ResultPath", "Output"))
            / fs::filesystem::path(solFilename);
    }

    env->report->outputOptionsReport();
    env->report->outputProblemInstanceReport();

    if(!solver.solveProblem()) // Solve the problem
    {
        env->output->outputCritical(" Error when solving problem.");
        return (0);
    }

    env->report->outputSolutionReport();

    env->output->outputInfo("╶──────────────────────────────────────────────────────────────────────────────────"
                            "───────────────────────────────────╴\r\n");

    std::string osrl = solver.getResultsOSrL();

    if(resultFile.empty())
    {
        fs::filesystem::path resultPath(env->settings->getSetting<std::string>("ResultPath", "Output"));
        resultPath /= env->settings->getSetting<std::string>("ProblemName", "Input");
        resultPath = resultPath.replace_extension(".osrl");

        if(!Utilities::writeStringToFile(resultPath.string(), osrl))
            env->output->outputCritical(" Error when writing OSrL file to: " + resultPath.string());
        else
            env->output->outputInfo(" Results written to: " + resultPath.string());
    }
    else
    {
        if(!Utilities::writeStringToFile(resultFile.string(), osrl))
            env->output->outputCritical(" Error when writing OSrL file to: " + resultFile.string());
        else
            env->output->outputInfo(" Results written to: " + resultFile.string());
    }

    if(cmdl["--trc"] || cmdl("--trc"))
    {
        std::string trace = solver.getResultsTrace();

        if(traceFile.empty())
        {
            fs::filesystem::path tracePath(env->settings->getSetting<std::string>("ResultPath", "Output"));
            tracePath /= env->settings->getSetting<std::string>("ProblemName", "Input");
            tracePath = tracePath.replace_extension(".trc");

            if(!Utilities::writeStringToFile(tracePath.string(), trace))
                env->output->outputCritical(" Error when writing trace file: " + tracePath.string());
            else
                env->output->outputInfo("                     " + tracePath.string());
        }
        else
        {
            if(!Utilities::writeStringToFile(traceFile.string(), trace))
                env->output->outputCritical(" Error when writing trace file: " + traceFile.string());
            else
                env->output->outputInfo("                     " + traceFile.string());
        }
    }

    if(cmdl["--sol"] || cmdl("--sol") || useASL)
    {
        std::string sol = solver.getResultsSol();

        if(solFile.empty())
        {
            fs::filesystem::path solPath(filename);
            solPath = solPath.replace_extension(".sol");

            if(!Utilities::writeStringToFile(solPath.string(), sol))
                env->output->outputCritical(" Error when writing AMPL sol file: " + solPath.string());
            else
                env->output->outputInfo("                     " + solPath.string());
        }
        else
        {
            if(!Utilities::writeStringToFile(solFile.string(), sol))
                env->output->outputCritical(" Error when writing AMPL sol file: " + solFile.string());
            else
                env->output->outputInfo("                     " + solFile.string());
        }
    }

    env->output->outputInfo("\r\n");
    env->output->outputInfo(" Log written to:     " + logFile.string());

    if(env->settings->getSetting<bool>("Debug.Enable", "Output"))
    {
        auto debugDirectory = fs::filesystem::current_path()
            / fs::filesystem::path(env->settings->getSetting<std::string>("Debug.Path", "Output"));
        env->output->outputInfo(" Debug directory:    " + debugDirectory.string());
    }

    env->results = NULL;

    return (0);
}