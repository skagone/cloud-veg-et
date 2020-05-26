import pkgutil

# this is the package we are inspecting -- for example 'email' from stdlib
import vegetLib

package = vegetLib
for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):
    print ("Found submodule %s " % (modname))
