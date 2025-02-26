#include "stl_pluginentry.h"
#include "stl_steppluginmanager.h"

ST_STL_PluginEntry::ST_STL_PluginEntry()
{
    _stepPluginManager = new ST_STL_StepPluginManager();
}

ST_STL_PluginEntry::~ST_STL_PluginEntry()
{
    delete _stepPluginManager;
}

QString ST_STL_PluginEntry::getVersion() const
{
    return "0.1";
}

CT_AbstractStepPlugin* ST_STL_PluginEntry::getPlugin() const
{
    return _stepPluginManager;
}

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
Q_EXPORT_PLUGIN2(plug_steplight, ST_STL_PluginEntry)
#endif
