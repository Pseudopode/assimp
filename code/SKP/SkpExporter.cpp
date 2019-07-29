#ifndef ASSIMP_BUILD_NO_EXPORT
#ifndef ASSIMP_BUILD_NO_SKP_EXPORTER

#include "SkpExporter.h"
#include <assimp/Bitmap.h>
#include <assimp/fast_atof.h>
#include <assimp/SceneCombiner.h>
#include <assimp/StringUtils.h>
#include <assimp/DefaultIOSystem.h>
#include <assimp/IOSystem.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>

#include <assimp/Exceptional.h>

#include <SketchUpAPI/sketchup.h>

#include <cassert>
#include <ctime>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

using namespace Assimp;

namespace Assimp {
namespace {

#define SU(api_function_call) {\
SUResult su_api_result = api_function_call;\
assert(SU_ERROR_NONE == su_api_result);\
}\

std::string GetString(const SUStringRef& string) {
    size_t length = 0;
    SU(SUStringGetUTF8Length(string, &length));
    std::vector<char> buffer(length + 1);
    size_t out_length = 0;
    SU(SUStringGetUTF8(string, length + 1, buffer.data(), &out_length));
    assert(out_length == length);
    return std::string(begin(buffer), end(buffer));
}

const double MeterToInch = 39.3700787402;


SUMaterialRef LoadMaterial(aiMaterial* ai_material, SUModelRef model) {
    SUMaterialRef material = SU_INVALID;
    SU(SUMaterialCreate(&material));

    std::cout << "Material: " << ai_material->GetName().C_Str() << "\n";
    SU(SUMaterialSetName(material, ai_material->GetName().C_Str())); // TODO: Ensure is unique?

    aiColor4D ai_color;
    if (aiReturn_SUCCESS == ai_material->Get(AI_MATKEY_COLOR_DIFFUSE, ai_color)) {
        SUByte r = static_cast<SUByte>(std::round(255.0 * ai_color.r));
        SUByte g = static_cast<SUByte>(std::round(255.0 * ai_color.g));
        SUByte b = static_cast<SUByte>(std::round(255.0 * ai_color.b));
        SUByte a = static_cast<SUByte>(std::round(255.0 * ai_color.a));
        std::cout << "  aiColor: " << ai_color.r << ", " << ai_color.g << ", " << ai_color.b << ", " << ai_color.a << "\n";
        std::cout << "  SUColor: " << unsigned(r) << ", " << unsigned(g) << ", " << unsigned(b) << ", " << unsigned(a) << "\n";

        SUColor color{ r, g, b, a };
        SU(SUMaterialSetColor(material, &color));
    }

    ai_real opacity;
    if (aiReturn_SUCCESS == ai_material->Get(AI_MATKEY_OPACITY, opacity)) {
        // TODO: SUMaterialSetUseOpacity ?
        SU(SUMaterialSetOpacity(material, opacity));
    }

    if (ai_material->GetTextureCount(aiTextureType_DIFFUSE) > 0) {
        aiString ai_path;
        aiTextureMapping ai_mapping; // Handle only aiTextureMapping_UV
        aiTextureMapMode ai_mode; // Handle only aiTextureMapMode_Wrap
        auto result = ai_material->GetTexture(aiTextureType_DIFFUSE, 0,
            &ai_path, &ai_mapping, nullptr, nullptr, nullptr, &ai_mode);
        assert(result == aiReturn_SUCCESS);
        assert(ai_mapping == aiTextureMapping_UV);
        assert(ai_mode == aiTextureMapMode_Wrap);
        std::cout << "  Texture: " << ai_path.C_Str() << "\n";

        double scale_s = 1.0;
        double scale_t = 1.0;
        SUTextureRef texture = SU_INVALID;
        SU(SUTextureCreateFromFile(&texture, ai_path.C_Str(), scale_s, scale_t)); // TODO: Handle missing file.
        SU(SUMaterialSetTexture(material, texture));
    }

    //SU(SUModelAddMaterials(model, 1, &material));
    auto result = SUModelAddMaterials(model, 1, &material);
    assert(result == SU_ERROR_NONE);
    return material;
}

void SetModelOptions(SUModelRef model) {
    SUOptionsManagerRef options = SU_INVALID;
    SU(SUModelGetOptionsManager(model, &options));

    SUOptionsProviderRef unit_options = SU_INVALID;
    SU(SUOptionsManagerGetOptionsProviderByName(options,
        "UnitsOptions", &unit_options));

    SUTypedValueRef value = SU_INVALID;
    SU(SUTypedValueCreate(&value));

    SU(SUTypedValueSetInt32(value, 0)); // 0: Decimal
    SU(SUOptionsProviderSetValue(unit_options, "LengthFormat", value));

    SU(SUTypedValueSetInt32(value, 2)); // 2: mm, 4: m
    SU(SUOptionsProviderSetValue(unit_options, "LengthUnit", value));

    SU(SUTypedValueSetInt32(value, 1)); // 1: 0.0
    SU(SUOptionsProviderSetValue(unit_options, "LengthPrecision", value));

    SU(SUTypedValueRelease(&value));
}

SUTransformation ColladaToSkpCoordinates() {
    SUTransformation scale_transformation;
    SU(SUTransformationScale(&scale_transformation, MeterToInch));

    SUPoint3D origin{ 0.0, 0.0, 0.0 };
    SUVector3D x_axis{ -1.0, 0.0, 0.0 }; // TODO: Is negative X-Axis correct?
    SUVector3D y_axis{ 0.0, 0.0, 1.0 };
    SUVector3D z_axis{ 0.0, 1.0, 0.0 };
    SUTransformation axes_transformation;
    SU(SUTransformationSetFromPointAndAxes(&axes_transformation,
        &origin, &x_axis, &y_axis, &z_axis));

    SUTransformation transformation;
    SU(SUTransformationMultiply(&axes_transformation, &scale_transformation,
        &transformation));

    return transformation;
}


class SkpExporter {
public:
    SkpExporter();

    void Export(const char* filepath, const aiScene* scene);

private:
    void LoadMaterials(const aiScene* scene);
    void LoadMeshDefinitions(const aiScene* scene);
    void LoadNodes(const aiScene* scene);

    SUComponentDefinitionRef MeshToDefinition(aiMesh* mesh);
    void NodeToInstance(aiNode* node, SUEntitiesRef entities);

    // Collada to SKP coordinate system transformation.
    const SUTransformation transformation_{};

    SUModelRef model_ = SU_INVALID;
    std::vector<SUMaterialRef> materials_;
    std::vector<SUComponentDefinitionRef> definitions_;
};


SkpExporter::SkpExporter() : transformation_(ColladaToSkpCoordinates()) {
}

void SkpExporter::Export(const char* filepath, const aiScene* scene) {
    assert(SUIsInvalid(model_));
    model_ = SU_INVALID;
    SU(SUModelCreate(&model_));
    SetModelOptions(model_);
    LoadMaterials(scene);
    LoadMeshDefinitions(scene);
    LoadNodes(scene);
    SU(SUModelMergeCoplanarFaces(model_)); // TODO: This might ignore materials...
    SU(SUModelSaveToFile(model_, filepath));
    SU(SUModelRelease(&model_));
}

void SkpExporter::LoadMaterials(const aiScene* scene) {
    for (size_t i = 0; i < scene->mNumMaterials; i++) {
        auto ai_material = scene->mMaterials[i];
        auto material = LoadMaterial(ai_material, model_);
        materials_.emplace_back(material);
    }
}

void SkpExporter::LoadMeshDefinitions(const aiScene* scene) {
    for (size_t i = 0; i < scene->mNumMeshes; i++) {
        auto mesh = scene->mMeshes[i];
        auto definition = MeshToDefinition(mesh);
        definitions_.emplace_back(definition);
    }
}

void SkpExporter::LoadNodes(const aiScene* scene) {
    SUEntitiesRef entities = SU_INVALID;
    SU(SUModelGetEntities(model_, &entities));
    NodeToInstance(scene->mRootNode, entities);
}

SUVector3D ComputeNormal(const aiMesh* mesh, const aiFace& face) {
    // Copied from GenFaceNormalsProcess.
    const aiVector3D* pV1 = &mesh->mVertices[face.mIndices[0]];
    const aiVector3D* pV2 = &mesh->mVertices[face.mIndices[1]];
    const aiVector3D* pV3 = &mesh->mVertices[face.mIndices[face.mNumIndices - 1]];
    const aiVector3D vNor = ((*pV2 - *pV1) ^ (*pV3 - *pV1)).NormalizeSafe();
    return { vNor.x, vNor.y, vNor.z };
}

SUComponentDefinitionRef SkpExporter::MeshToDefinition(aiMesh* mesh) {
    assert(mesh->HasFaces());
    assert(mesh->HasPositions());
    assert(mesh->HasNormals());

    //std::cout << "Mesh: " << mesh->mName.C_Str() << "\n";

    SUMaterialRef material = materials_.at(mesh->mMaterialIndex);

    SUGeometryInputRef input = SU_INVALID;
    SU(SUGeometryInputCreate(&input));

    for (size_t i = 0; i < mesh->mNumVertices; i++) {
        auto vertex = mesh->mVertices[i];
        SUPoint3D position{ vertex.x, vertex.y, vertex.z };
        SU(SUPoint3DTransform(&transformation_, &position));
        SU(SUGeometryInputAddVertex(input, &position));
    }

    for (size_t i = 0; i < mesh->mNumFaces; i++) {
        auto face = mesh->mFaces[i];
        //std::cout << "  Face: " << i << "\n";

        SULoopInputRef loop = SU_INVALID;
        SU(SULoopInputCreate(&loop));

        for (size_t j = 0; j < face.mNumIndices; j++) {
            size_t vertex_index = face.mIndices[j];
            SU(SULoopInputAddVertexIndex(loop, vertex_index));
        }

        // Set edges to soft+smooth if vertex normals isn't in the same
        // direction as the face normal. This is the best SketchUp can do since
        // it doesn't allow explicit control over vertex normals.
        SUVector3D face_normal = ComputeNormal(mesh, face);
        //std::cout << "    Normal: " << face_normal.x << ", " << face_normal.y << ", " << face_normal.z << "\n";
        for (size_t j = 0; j < face.mNumIndices; j++) {
            // Get the index of the next vertex in the loop.
            size_t k = (j + 1) % face.mNumIndices;
            // Fetch the vertex normals.
            size_t ni1 = face.mIndices[j];
            size_t ni2 = face.mIndices[k];
            auto n1 = mesh->mNormals[ni1];
            auto n2 = mesh->mNormals[ni2];
            SUVector3D normal1{ n1.x, n1.y, n1.z };
            SUVector3D normal2{ n2.x, n2.y, n2.z };
            // Check if they are in the same direction as the face.
            bool same1 = false, same2 = false;
            SU(SUVector3DIsSameDirectionAs(&face_normal, &normal1, &same1));
            SU(SUVector3DIsSameDirectionAs(&face_normal, &normal2, &same2));

            //std::cout << "      V1 Normal: " << normal1.x << ", " << normal1.y << ", " << normal1.z << " same: " << same1 << "\n";
            //std::cout << "      V2 Normal: " << normal2.x << ", " << normal2.y << ", " << normal2.z << " same: " << same2 << "\n";

            // If either vertex normal is not in the same direction then assume
            // this edge should be soft+smooth.
            if (!same1 || !same2) {
                //std::cout << "    Soft+Smooth edge: " << j << "\n";
                size_t edge_index = j;
                SU(SULoopInputEdgeSetSoft(loop, edge_index, true));
                // TODO: SketchUp Bug! The smooth propery appear to be applied
                //   to the whole loop instead of just the given edge index.
                //SU(SULoopInputEdgeSetSmooth(loop, edge_index, true));
            }
        }

        size_t face_index = 0;
        SU(SUGeometryInputAddFace(input, &loop, &face_index));

        size_t num_uv_coords = 0;
        SUMaterialInput material_input{ num_uv_coords, {}, {}, material };
        SU(SUGeometryInputFaceSetFrontMaterial(input, face_index, &material_input));
    }

    SUComponentDefinitionRef definition = SU_INVALID;
    SU(SUComponentDefinitionCreate(&definition));

    auto name = mesh->mName;
    SU(SUComponentDefinitionSetName(definition, name.C_Str()));

    // Add the definition to the model before populating it with geometry.
    // Otherwise it'll lead PIDs not being generated and yield warnings when
    // the user loads the model. When materials have been assigned to entities
    // in the definition it all lead to random run-time crashes.
    SU(SUModelAddComponentDefinitions(model_, 1, &definition));

    SUEntitiesRef entities = SU_INVALID;
    SU(SUComponentDefinitionGetEntities(definition, &entities));

    SU(SUEntitiesFill(entities, input, true));
    SU(SUGeometryInputRelease(&input));

    // Fix the smooth bug in SULoopInputEdgeSetSmooth.
    size_t num_edges = 0;
    SU(SUEntitiesGetNumEdges(entities, false, &num_edges));
    std::vector<SUEdgeRef> edges(num_edges, SU_INVALID);
    SU(SUEntitiesGetEdges(entities, false, num_edges, edges.data(), &num_edges));
    for (auto& edge : edges) {
        bool is_soft = false;
        SU(SUEdgeGetSoft(edge, &is_soft));
        if (is_soft) {
            SU(SUEdgeSetSmooth(edge, true));
        }
    }

    return definition;
}

void SkpExporter::NodeToInstance(aiNode* node, SUEntitiesRef entities) {
    std::cout << "Node: " << node->mName.C_Str() << "\n";

    if (node->mNumChildren == 0 && node->mNumMeshes == 0) {
        return;
    }

    auto tr = node->mTransformation;
    std::cout << "  Transformation:\n";
    std::cout << tr.a1 << ", " << tr.a2 << ", " << tr.a3 << ", " << tr.a4 << "\n";
    std::cout << tr.b1 << ", " << tr.b2 << ", " << tr.b3 << ", " << tr.b4 << "\n";
    std::cout << tr.c1 << ", " << tr.c2 << ", " << tr.c3 << ", " << tr.c4 << "\n";
    std::cout << tr.d1 << ", " << tr.d2 << ", " << tr.d3 << ", " << tr.d4 << "\n";

    // TODO: Verify
    // assimp is row-major while SketchUp is column-major
    //const double scale = MeterToInch;
    //SUTransformation transformation{
    //    tr.a1, tr.b1, tr.c1, tr.d1,
    //    tr.a2, tr.b2, tr.c2, tr.d2,
    //    tr.a3, tr.b3, tr.c3, tr.d3,
    //    tr.a4 * scale, tr.b4 * scale, tr.c4 * scale, tr.d4,
    //};

    const double scale = MeterToInch;
    SUTransformation transformation{
        tr.a1, tr.b1, tr.c1, tr.d1,
        tr.a2, tr.b2, tr.c2, tr.d2,
        tr.a3, tr.b3, tr.c3, tr.d3,
        // Beware that the Y and Z position is swapped to account for the YZ
        // axis difference.
        tr.a4 * scale, tr.c4 * scale, tr.b4 * scale, tr.d4,
    };

    std::cout << "  Meshes: " << node->mNumMeshes << "\n";
    for (size_t i = 0; i < node->mNumMeshes; i++) {
        auto mesh_index = node->mMeshes[i];
        auto definition = definitions_.at(mesh_index);
        std::cout << "    MeshIndex: " << mesh_index << "\n";

        SUComponentInstanceRef instance = SU_INVALID;
        SU(SUComponentDefinitionCreateInstance(definition, &instance));
        SU(SUComponentInstanceSetName(instance, node->mName.C_Str()));
        SU(SUComponentInstanceSetTransform(instance, &transformation));

        SU(SUEntitiesAddInstance(entities, instance, nullptr));
    }

    for (size_t i = 0; i < node->mNumChildren; i++) {
        auto child = node->mChildren[i];

        SUGroupRef group = SU_INVALID;
        SU(SUGroupCreate(&group));
        SU(SUGroupSetName(group, child->mName.C_Str()));
        SU(SUEntitiesAddGroup(entities, group));

        SUEntitiesRef child_entities = SU_INVALID;
        SU(SUGroupGetEntities(group, &child_entities));

        NodeToInstance(child, child_entities);
    }
}

} // namespace


// ------------------------------------------------------------------------------------------------
// Worker function for exporting a scene to Collada. Prototyped and registered in Exporter.cpp
void ExportSceneSKP(const char* pFile, IOSystem* pIOSystem,
        const aiScene* pScene, const ExportProperties* /*pProperties*/) {

    // <debug>
    std::cout << "\n";
    std::cout << "Metadata:\n";
    for (size_t i = 0; i < pScene->mMetaData->mNumProperties; i++) {
        auto key = pScene->mMetaData->mKeys[i];
        //auto value = pScene->mMetaData->mValues[i];
        std::cout << key.C_Str() << "\n";
    }
    std::cout << "\n";
    // </debug>

    SUInitialize();
    {
        SkpExporter exporter;
        exporter.Export(pFile, pScene);
    }
    SUTerminate();
}

} // end of namespace Assimp

#endif
#endif
