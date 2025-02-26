#ifndef ST_STL_PLUGINENTRY_H
#define ST_STL_PLUGINENTRY_H

#include "pluginentryinterface.h"

class ST_STL_StepPluginManager;

class ST_STL_PluginEntry : public PluginEntryInterface
{
    Q_OBJECT

#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
    Q_PLUGIN_METADATA(IID PluginEntryInterface_iid)
#endif

    Q_INTERFACES(PluginEntryInterface)

public:
    ST_STL_PluginEntry();
    ~ST_STL_PluginEntry() override;

    QString getVersion() const override;
    CT_AbstractStepPlugin* getPlugin() const override;

private:
    ST_STL_StepPluginManager *_stepPluginManager;
};


#endif // ST_STL_PLUGINENTRY_H
